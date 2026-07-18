#!/usr/bin/env python3
"""Hash-checked loading boundary for a future physical SU2R CMB backend."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Mapping

import numpy as np

from cmb_theory_contract import (
    BackgroundPrediction,
    CmbSpectra,
    CmbTheoryAdapter,
    FloatArray,
    SpectralDistortion,
    TheoryCapabilities,
)


API_VERSION = 1
LENSING_EQUATIONS = (
    "action_or_field_equations",
    "background_field_equations",
    "scalar_dynamical_equations",
    "scalar_constraint_equations",
    "perturbed_stress_energy",
    "gauge_and_variable_dictionary",
    "initial_conditions",
    "weyl_potential_and_lensing_relation",
)
PRIMARY_CMB_EQUATIONS = LENSING_EQUATIONS + (
    "photon_baryon_neutrino_couplings",
    "recombination_and_thermal_history",
)


class PhysicalTheoryIncomplete(RuntimeError):
    """Raised when a physical-model label would exceed the supplied equations."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class SU2RRegistry:
    path: Path
    payload: dict[str, object]

    @classmethod
    def load(cls, path: Path) -> "SU2RRegistry":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if int(payload.get("api_version", -1)) != API_VERSION:
            raise ValueError(f"SU2R registry must use API version {API_VERSION}")
        return cls(path=path.resolve(), payload=payload)

    def missing(self, branch: str) -> list[str]:
        if branch == "lensing":
            required = LENSING_EQUATIONS
        elif branch == "primary_cmb":
            required = PRIMARY_CMB_EQUATIONS
        else:
            raise ValueError(f"Unknown SU2R branch: {branch}")
        equations = self.payload.get("equations", {})
        if not isinstance(equations, dict):
            return list(required)
        missing = [name for name in required if not equations.get(name)]
        implementation = self.payload.get("implementation", {})
        if not isinstance(implementation, dict):
            return missing + ["implementation.module_path", "implementation.sha256"]
        if not implementation.get("module_path"):
            missing.append("implementation.module_path")
        if not implementation.get("sha256"):
            missing.append("implementation.sha256")
        return missing

    def require(self, branch: str) -> None:
        missing = self.missing(branch)
        if missing:
            raise PhysicalTheoryIncomplete(
                f"Physical SU2R {branch} backend is incomplete: {', '.join(missing)}"
            )


def _load_backend(registry: SU2RRegistry) -> ModuleType:
    implementation = registry.payload["implementation"]
    assert isinstance(implementation, dict)
    raw_path = Path(str(implementation["module_path"]))
    module_path = raw_path if raw_path.is_absolute() else registry.path.parent / raw_path
    module_path = module_path.resolve()
    expected_hash = str(implementation["sha256"]).lower()
    measured_hash = sha256_file(module_path)
    if measured_hash != expected_hash:
        raise RuntimeError(
            f"SU2R backend hash mismatch: expected {expected_hash}, measured {measured_hash}"
        )
    spec = importlib.util.spec_from_file_location("registered_su2r_physical_backend", module_path)
    if spec is None or spec.loader is None:
        raise ImportError(module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    for name in ("background", "spectra"):
        if not callable(getattr(module, name, None)):
            raise TypeError(f"Registered SU2R backend lacks callable {name}()")
    return module


class RegisteredSU2RPhysicalAdapter(CmbTheoryAdapter):
    """Delegate to a backend only after equation and source-hash validation."""

    def __init__(self, registry_path: Path, branch: str = "lensing") -> None:
        self.registry = SU2RRegistry.load(registry_path)
        self.registry.require(branch)
        self.branch = branch
        self.backend = _load_backend(self.registry)
        self.capabilities = TheoryCapabilities(
            background=True,
            thermodynamics=branch == "primary_cmb",
            scalar_perturbations=True,
        )

    def background(self, parameters: Mapping[str, float], z: FloatArray) -> BackgroundPrediction:
        prediction = self.backend.background(parameters, np.asarray(z, dtype=float))
        if not isinstance(prediction, BackgroundPrediction):
            raise TypeError("SU2R backend background() returned the wrong contract type")
        prediction.validate()
        return prediction

    def spectra(self, parameters: Mapping[str, float], lmax: int) -> CmbSpectra:
        spectra = self.backend.spectra(parameters, int(lmax))
        if not isinstance(spectra, CmbSpectra):
            raise TypeError("SU2R backend spectra() returned the wrong contract type")
        spectra.validate(("TT", "TE", "EE", "BB", "PP"))
        return spectra

    def distortion(self, parameters: Mapping[str, float], frequency_ghz: FloatArray) -> SpectralDistortion:
        function = getattr(self.backend, "distortion", None)
        if not callable(function):
            return super().distortion(parameters, frequency_ghz)
        prediction = function(parameters, np.asarray(frequency_ghz, dtype=float))
        if not isinstance(prediction, SpectralDistortion):
            raise TypeError("SU2R backend distortion() returned the wrong contract type")
        prediction.validate()
        return prediction
