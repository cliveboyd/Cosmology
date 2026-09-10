from __future__ import annotations

import csv
import json
import math
import shutil
from datetime import date
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


DOCS_ROOT = Path(r"C:\Users\clive\Documents\Codex")
COSMO_ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
OUTPUT_ROOT = Path(
    r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs"
    r"\bbn_lithium_full_cmb_nuisance_20260821"
)
POST_ROOT = OUTPUT_ROOT / "postprocess"
QUAL_ROOT = OUTPUT_ROOT / "qualification_audit_combined_ess_topup"
GITHUB_PAPERS = COSMO_ROOT / "github_export" / "docs" / "papers"
GITHUB_CODE = COSMO_ROOT / "github_export" / "code" / "bbn_lithium"

REPORT_PATH = DOCS_ROOT / "big_bang_nucleosynthesis_hydrogen_helium_lithium_report_2026-09-11.docx"
PAPER_PATH = DOCS_ROOT / "bbn_lithium_august_findings_formal_scientific_paper_2026-09-11.docx"

PROJECT_URL = "https://github.com/cliveboyd/Cosmology"
BBN_CODE_URL = PROJECT_URL + "/tree/main/github_export/code/bbn_lithium"
BBN_PAPERS_URL = PROJECT_URL + "/tree/main/github_export/docs/papers"


REFERENCES = [
    "Particle Data Group, Big-Bang Nucleosynthesis review, Review of Particle Physics (2024).",
    "B. D. Fields, The Primordial Lithium Problem, Annual Review of Nuclear and Particle Science 61, 47 (2011).",
    "R. H. Cyburt, B. D. Fields and K. A. Olive, A Bitter Pill: The Primordial Lithium Problem Worsens, JCAP 11, 012 (2008), arXiv:0808.2818.",
    "R. H. Cyburt et al., Big Bang Nucleosynthesis: 2015, Rev. Mod. Phys. 88, 015004.",
    "B. D. Fields et al., Big-Bang Nucleosynthesis After Planck, JCAP 03 (2020) 010.",
    "C. Pitrou, A. Coc, J.-P. Uzan and E. Vangioni, Precision Big Bang Nucleosynthesis with Improved Helium-4 Predictions, Phys. Rep. 754, 1-66 (2018).",
    "F. Spite and M. Spite, Abundance of lithium in unevolved halo stars and old disc stars, A&A 115, 357 (1982).",
    "M. Spite and F. Spite, Lithium abundance at the formation of the Galaxy, Nature 297, 483-485 (1982).",
    "L. Sbordone et al., The metal-poor end of the Spite plateau. I. Stellar parameters, metallicities, and lithium abundances, A&A 522, A26 (2010).",
    "S. G. Ryan, A critique of the Spite Plateau and the astration of primordial lithium, MNRAS 522, 1358-1376 (2023).",
    "R. J. Cooke et al., One percent determination of the primordial deuterium abundance, ApJ 855, 102.",
    "M. L. Giannotti et al., LINX: A Fast, Differentiable and Extensible Big Bang Nucleosynthesis Package (2024).",
    "Joint LINX collaboration, Cosmological Parameter Estimation with a Joint-Likelihood Analysis of the CMB and BBN (2024).",
    "P. F. Depta et al., ACROPOLIS: A generic Framework for Photodisintegration of Light elements, JCAP 03 (2021) 061.",
    "P. F. Depta et al., Updated BBN constraints on electromagnetic decays of MeV-scale particles (2020).",
    "K. Langhoff et al., New insights into axion freeze-in, JHEP 11 (2024) 166.",
    "E. Braaten and D. Segel, Neutrino energy loss from the plasma process at all temperatures and densities, Phys. Rev. D 48, 1478.",
    "D. J. Fixsen et al., The Cosmic Microwave Background Spectrum from the Full COBE FIRAS Data Set, ApJ 473, 576.",
    "R. Khatri and R. A. Sunyaev, Creation of the CMB spectrum: precise analytic solutions for the blackbody photosphere, JCAP 06 (2012) 038.",
    "Planck Collaboration, Planck 2018 results. VI. Cosmological parameters, A&A 641, A6 (2020).",
    "G. Efstathiou and S. Gratton, A detailed description of the CamSpec likelihood pipeline and a reanalysis of the Planck high frequency maps, MNRAS 496, L91 (2020); PR4 update.",
    "J. Carron et al., Planck 2018 CMB lensing reconstruction from the final NPIPE data release.",
    "Spectroxide project implementation, numerical thermalisation and spectral-distortion framework (2026).",
]


ACRONYMS = [
    ("ACROPOLIS", "Public electromagnetic-cascade framework used for BBN photodisintegration constraints."),
    ("AIC", "Akaike Information Criterion; lower is preferred after a penalty of 2 parameters."),
    ("BBN", "Big Bang Nucleosynthesis."),
    ("BIC", "Bayesian Information Criterion; lower is preferred with a stronger log(N) parameter penalty."),
    ("CAMB", "Code for Anisotropies in the Microwave Background."),
    ("CMB", "Cosmic microwave background."),
    ("Cobaya", "Cosmological Bayesian analysis package used for the Planck nuisance MCMC."),
    ("D/H", "Deuterium-to-hydrogen abundance ratio."),
    ("ESS", "Effective sample size; here the autocorrelation-aware number of effectively independent posterior samples."),
    ("FIRAS", "Far Infrared Absolute Spectrophotometer on COBE."),
    ("LINX", "Fast differentiable BBN package used for abundance calculations and validation."),
    ("MCMC", "Markov-chain Monte Carlo."),
    ("NPIPE", "Planck reprocessing pipeline used for PR4 products."),
    ("PR4", "Planck public release 4."),
    ("SPD", "Symmetric positive definite; used for covariance handling and Cholesky/LAPACK inversion."),
    ("TV index", "Traceability/validation index. The reports place it after each investigation section."),
]


VARIABLES = [
    ("eta10", "10^10 times the baryon-to-photon ratio eta; converted from omega_b using eta10 = 273.78 omega_b for T_CMB = 2.7255 K."),
    ("omega_b / ombh2", "Physical baryon density Omega_b h^2."),
    ("omega_c / omch2", "Physical cold dark matter density Omega_c h^2."),
    ("N_eff / nnu", "Effective relativistic species count inferred independently in the CMB nuisance chain."),
    ("Y_p / YHe", "Primordial helium-4 mass fraction; sampled independently in the full CMB likelihood."),
    ("D/H", "Deuterium abundance anchor; required to remain within the adopted observational gate."),
    ("Li-7/H", "Lithium-7 abundance. Most primordial mass-7 is produced as Be-7 and later decays to Li-7."),
    ("m_phi", "Mediator mass in MeV for the selective electromagnetic-transfer candidate."),
    ("E_gamma", "Characteristic photon energy, approximated as m_phi/2 for two-photon decay."),
    ("tau_phi", "Mediator lifetime in seconds."),
    ("n_phi/n_gamma", "Mediator abundance relative to photons at T = 10 MeV."),
    ("g_phi_gamma", "Effective mediator-photon coupling in GeV^-1."),
    ("mu", "CMB spectral-distortion chemical-potential style readout."),
    ("R-1", "Convergence diagnostic used by Cobaya and the split-chain audit."),
    ("chi2", "Goodness-of-fit statistic; lower values indicate a better fit for the same data and model complexity."),
]


PROGRAMS = [
    ("analyze_bbn_lithium_linx_fr_network_2026-07-16.py", "July LINX/FR network controls; reproduced the lithium problem and did not identify a D+He+Li pass."),
    ("analyze_su2_bbn_lithium_gate_2026-07-17.py", "July SU2-style expansion/thermal gate for lithium; part of the negative control phase."),
    ("validate_bbn_lithium_linx_matrix_2026-07-17.py", "Strict July validation matrix for the BBN lithium controls."),
    ("audit_bbn_lithium_linx_fr_scan_2026-07-17.py", "Audit layer for the July LINX/FR scan."),
    ("run_cobaya_spd_cache_builder_2026_08_21.py", "Built the validated SPD cache for the Planck PR4/NPIPE CamSpec likelihood."),
    ("postprocess_planck_pr4_full_nuisance.py", "Post-processed the native full CMB nuisance chain and compared the registered lithium cell with the CMB posterior."),
    ("audit_planck_pr4_full_nuisance_qualification.py", "Strict qualification audit: split chains, autocorrelation-aware ESS, registered-cell density, nuisance stability and half-chain stability."),
    ("start_planck_pr4_full_nuisance_background.ps1", "Background launcher for the native Planck PR4 full nuisance MCMC."),
    ("start_planck_pr4_full_nuisance_ess_topup_background.ps1", "Background launcher for the ESS top-up chain."),
    ("build_bbn_li_reports_20260911.py", "This reproducible DOCX report generator."),
]


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def fmt(value, digits: int = 4) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        try:
            value = float(value)
        except ValueError:
            return value
    if isinstance(value, bool):
        return "PASS" if value else "FAIL"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        if value == 0:
            return "0"
        if abs(value) >= 1e5 or abs(value) < 1e-3:
            return f"{value:.{digits}e}"
        return f"{value:.{digits}f}".rstrip("0").rstrip(".")
    return str(value)


def set_cell_fill(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), fill)
    tc_pr.append(shd)


def set_cell_text(cell, text: str, bold: bool = False) -> None:
    cell.text = ""
    p = cell.paragraphs[0]
    r = p.add_run(str(text))
    r.bold = bold
    r.font.size = Pt(8.5)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.TOP


def add_table(doc: Document, headers: list[str], rows: list[list[str]], title: str | None = None) -> None:
    if title:
        p = doc.add_paragraph()
        p.style = "Caption"
        p.add_run(title).bold = True
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, h in enumerate(headers):
        set_cell_text(table.rows[0].cells[i], h, True)
        set_cell_fill(table.rows[0].cells[i], "D9EAF7")
    for row in rows:
        cells = table.add_row().cells
        for i, val in enumerate(row):
            set_cell_text(cells[i], val)
    doc.add_paragraph()


def add_title_block(doc: Document, title: str, subtitle: str, evidence_status: str) -> None:
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(title)
    r.bold = True
    r.font.size = Pt(18)
    r.font.color.rgb = RGBColor(31, 78, 121)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(subtitle)
    r.italic = True
    r.font.size = Pt(11)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run("Authors: Clive Stewart Boyd; David Ng").bold = True
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run("Generated: 11 September 2026 | Melbourne, Australia")
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run("Computational assistance and reproducibility support: ChatGPT Codex sol 5.5").italic = True
    doc.add_paragraph()
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("EVIDENCE STATUS: " + evidence_status)
    r.bold = True
    r.font.color.rgb = RGBColor(156, 87, 0)
    doc.add_page_break()


def add_para(doc: Document, text: str, style: str | None = None) -> None:
    p = doc.add_paragraph(text)
    if style:
        p.style = style


def add_bullets(doc: Document, items: list[str]) -> None:
    for item in items:
        p = doc.add_paragraph(style="List Bullet")
        p.add_run(item)


def add_numbered(doc: Document, items: list[str]) -> None:
    for item in items:
        p = doc.add_paragraph(style="List Number")
        p.add_run(item)


def add_tv_index(doc: Document, rows: list[tuple[str, str]]) -> None:
    add_table(doc, ["Traceability item", "Recorded value"], [[a, b] for a, b in rows], "TV index")


def setup_doc(title: str) -> Document:
    doc = Document()
    section = doc.sections[0]
    section.top_margin = Inches(0.55)
    section.bottom_margin = Inches(0.55)
    section.left_margin = Inches(0.62)
    section.right_margin = Inches(0.62)
    styles = doc.styles
    styles["Normal"].font.name = "Aptos"
    styles["Normal"].font.size = Pt(9.5)
    styles["Title"].font.name = "Aptos Display"
    styles["Heading 1"].font.name = "Aptos Display"
    styles["Heading 1"].font.size = Pt(15)
    styles["Heading 2"].font.name = "Aptos"
    styles["Heading 2"].font.size = Pt(12)
    styles["Heading 3"].font.name = "Aptos"
    styles["Heading 3"].font.size = Pt(10.5)
    for style_name in ["Caption", "List Bullet", "List Number"]:
        styles[style_name].font.name = "Aptos"
        styles[style_name].font.size = Pt(8.5 if style_name == "Caption" else 9.5)
    doc.core_properties.title = title
    doc.core_properties.author = "Clive Stewart Boyd; David Ng"
    return doc


def add_footer(doc: Document) -> None:
    for sec in doc.sections:
        footer = sec.footer.paragraphs[0]
        footer.text = "BBN-Li investigation | Clive Stewart Boyd; David Ng | Generated 11 September 2026"
        footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
        for run in footer.runs:
            run.font.size = Pt(8)


def gate_rows() -> list[list[str]]:
    rows = read_csv_dicts(QUAL_ROOT / "planck_pr4_full_nuisance_qualification_gates.csv")
    out = []
    for r in rows:
        out.append([r.get("gate", ""), r.get("threshold", ""), fmt(r.get("observed", "")), r.get("pass", "")])
    return out


def metric_rows() -> list[list[str]]:
    rows = read_csv_dicts(QUAL_ROOT / "planck_pr4_full_nuisance_audit_metrics.csv")
    keep = ["ombh2", "eta10", "nnu", "YHe", "A_planck", "amp_143", "amp_217", "amp_143x217", "calTE", "calEE"]
    out = []
    for r in rows:
        if r.get("parameter") in keep:
            out.append([
                r.get("parameter", ""),
                fmt(r.get("split_Rminus1", ""), 4),
                fmt(r.get("autocorr_ess", ""), 1),
                fmt(r.get("tau_int", ""), 1),
                fmt(r.get("first_last_mean_shift_sigma", ""), 3),
                r.get("status", ""),
            ])
    return out


def posterior_rows(selection: str = "all_rows") -> list[list[str]]:
    rows = read_csv_dicts(POST_ROOT / "planck_pr4_full_nuisance_posterior_summary.csv")
    keep = ["ombh2", "eta10", "nnu", "YHe", "H0", "omegam", "sigma8", "S8", "chi2__CMB"]
    out = []
    for r in rows:
        if r.get("selection") == selection and r.get("param") in keep:
            out.append([
                r.get("param", ""),
                fmt(r.get("mean", ""), 6),
                fmt(r.get("std", ""), 4),
                fmt(r.get("q50", ""), 6),
                f"{fmt(r.get('q16', ''), 6)} to {fmt(r.get('q84', ''), 6)}",
                f"{fmt(r.get('q2p5', ''), 6)} to {fmt(r.get('q97p5', ''), 6)}",
            ])
    return out


def candidate_rows() -> list[list[str]]:
    data = read_json(POST_ROOT / "planck_pr4_full_nuisance_bestfit.json")
    candidate = {
        "mediator_mass_MeV": 4.44,
        "photon_energy_MeV": 2.22,
        "lifetime_s": 158489.3,
        "n_phi_over_n_gamma_at_10MeV": 1.211528e-06,
        "g_phi_gamma_GeV_inv": 9.767274e-11,
        "D_H_times_1e5": 2.493161,
        "Yp": 0.24673503,
        "Li7_H_times_1e10": 1.488554,
        "He3_D": 0.415242,
        "eta10_CMB": 6.096517,
        "Neff_CMB": 3.044601,
        "mu_distortion": 1.476434e-09,
        "joint_pass_probability": 0.650158,
        "delta_log_evidence_vs_standard": 105.206012,
    }
    out = [[k, fmt(v, 6)] for k, v in candidate.items()]
    if data:
        out.append(["Planck bestfit eta10", fmt(data.get("eta10"), 6)])
        out.append(["Planck bestfit nnu", fmt(data.get("nnu"), 6)])
        out.append(["Planck bestfit YHe", fmt(data.get("YHe"), 6)])
        out.append(["Planck bestfit chi2__CMB", fmt(data.get("chi2__CMB"), 6)])
    return out


def registered_density_rows() -> list[list[str]]:
    d = read_json(QUAL_ROOT / "planck_pr4_full_nuisance_registered_density.json")
    return [
        ["posterior mean eta10", fmt(d.get("posterior_mean", {}).get("eta10"), 6)],
        ["posterior mean N_eff", fmt(d.get("posterior_mean", {}).get("nnu"), 6)],
        ["posterior mean Y_p", fmt(d.get("posterior_mean", {}).get("YHe"), 6)],
        ["registered eta10", fmt(d.get("registered_cell", {}).get("eta10"), 6)],
        ["registered N_eff", fmt(d.get("registered_cell", {}).get("nnu"), 6)],
        ["registered Y_p", fmt(d.get("registered_cell", {}).get("YHe"), 6)],
        ["registered Delta chi2", fmt(d.get("registered_delta_chi2_gaussian"), 6)],
        ["Mahalanobis distance", fmt(d.get("registered_mahalanobis_distance"), 6)],
        ["density ratio to posterior peak", fmt(d.get("registered_density_ratio_to_peak_gaussian"), 6)],
        ["posterior mass inside equal-density contour", fmt(d.get("posterior_mass_inside_registered_equal_density_contour"), 6)],
        ["joint 95 gate pass", str(d.get("joint_95_gate_pass"))],
        ["weighted kNN density ratio", fmt(d.get("weighted_knn_density_ratio_registered_to_peak"), 6)],
    ]


def add_image_if_exists(doc: Document, path: Path, caption: str, width: float = 6.7) -> None:
    if not path.exists():
        return
    doc.add_picture(str(path), width=Inches(width))
    p = doc.add_paragraph(caption)
    p.style = "Caption"


def add_common_appendices(doc: Document) -> None:
    doc.add_heading("Appendix A. Acronyms", level=1)
    add_table(doc, ["Acronym", "Meaning"], [[a, b] for a, b in ACRONYMS], "Table A1. Acronym glossary")
    doc.add_heading("Appendix B. Variables", level=1)
    add_table(doc, ["Variable", "Meaning"], [[a, b] for a, b in VARIABLES], "Table B1. Variable glossary")
    doc.add_heading("Appendix C. Program and Data Manifest", level=1)
    add_table(doc, ["Program or path", "Role in the investigation"], PROGRAMS, "Table C1. Principal programs")
    add_table(
        doc,
        ["Object", "Path or URL"],
        [
            ["Local output root", str(OUTPUT_ROOT)],
            ["Combined qualification audit", str(QUAL_ROOT)],
            ["Post-process folder", str(POST_ROOT)],
            ["GitHub project", PROJECT_URL],
            ["GitHub BBN code", BBN_CODE_URL],
            ["GitHub papers", BBN_PAPERS_URL],
        ],
        "Table C2. Repository and output cross-reference",
    )
    doc.add_heading("Appendix D. References", level=1)
    for idx, ref in enumerate(REFERENCES, start=1):
        doc.add_paragraph(f"{idx}. {ref}")


def add_chronology(doc: Document) -> None:
    doc.add_heading("Chronological Development and Test History", level=1)
    add_table(
        doc,
        ["Date", "Investigation stage", "Result status"],
        [
            [
                "16-17 July 2026",
                "Standard BBN, FR/clock-proxy, neutron-lifetime, reaction-rate and SU2-style expansion controls with LINX validation.",
                "Negative result: no joint D/H, helium and Li-7 solution. This is the explicit failure boundary for the preliminary July programme.",
            ],
            [
                "17-18 July 2026",
                "Hierarchical and surrogate follow-up tests, plus strict matrix validation.",
                "Useful for scoping but not a solution. Lithium could be moved only at the cost of other abundance or CMB gates.",
            ],
            [
                "20-21 August 2026",
                "Selective post-BBN electromagnetic-transfer hypothesis. Threshold screen, ACROPOLIS cascade tests, frozen LINX priors and production modelling.",
                "Positive candidate emerged: a narrow MeV-scale window can reduce mass-7 while preserving D/H, Y_p, He3/D, entropy and FIRAS gates.",
            ],
            [
                "21 August 2026",
                "Full Planck PR4/NPIPE CMB nuisance-likelihood protocol launched using Cobaya/CAMB and native Planck likelihoods.",
                "Initial status pending. The goal was to replace a compressed CMB Gaussian gate with a native spectrum-level posterior.",
            ],
            [
                "2 September 2026",
                "Main full CMB nuisance chain completed and post-processed.",
                "Positive but not yet fully qualified: Cobaya convergence was reported, but the stricter independent audit found autocorrelation ESS below the preregistered 1000 gate for some nuisance parameters.",
            ],
            [
                "2-11 September 2026",
                "ESS top-up chain run after the main chain to strengthen nuisance-tail sampling.",
                "Computer reset stopped the top-up, but not before adding enough samples for the combined audit.",
            ],
            [
                "11 September 2026",
                "Combined main plus ESS top-up qualification audit.",
                "Pass: all preregistered audit gates passed, including split-chain stability, autocorrelation-aware ESS, nuisance stability, half-chain stability and registered-cell posterior density.",
            ],
        ],
        "Table 1. Investigation chronology",
    )
    add_tv_index(
        doc,
        [
            ("Chronology source", "July 2026 local BBN/Li outputs, August full CMB nuisance folder, September combined qualification audit."),
            ("Failure boundary", "July preliminary programme failed to solve lithium while preserving D/H, helium and CMB-compatible baryon density."),
            ("Positive boundary", "August selective-transfer candidate became positive only after moving from global expansion/clock controls to a narrow post-BBN electromagnetic mechanism."),
        ],
    )


def build_update_report() -> None:
    doc = setup_doc("BBN hydrogen helium lithium report 2026-09-11")
    add_title_block(
        doc,
        "Normal Big Bang Nucleosynthesis and the Primordial Lithium Problem",
        "Hydrogen, helium and lithium production from standard-network controls to a qualified finite-temperature selective-transfer candidate",
        "Qualification-audit passed for the full Planck PR4/NPIPE CMB nuisance gate; scientific promotion remains conditional.",
    )

    doc.add_heading("Executive Summary", level=1)
    add_bullets(
        doc,
        [
            "Standard BBN remains successful for deuterium and helium at the adopted CMB-compatible baryon density, but predicts Li-7/H roughly three to four times above the warm metal-poor stellar plateau.",
            "The July 2026 investigations were negative: FR/clock-proxy changes, SU2-style expansion controls, reaction-rate pulls and network scans did not produce a joint D/H, helium and lithium solution.",
            "The August 2026 programme identified a selective post-BBN electromagnetic-transfer candidate. The candidate is narrow and mechanism-dependent, but it has passed the abundance, cascade, entropy, FIRAS and leading finite-temperature checks recorded so far.",
            "The new September 2026 result is important: the combined native Planck PR4/NPIPE nuisance audit passes every preregistered gate. The registered lithium cell is inside the joint 95% CMB posterior with Delta chi2 = 1.257.",
            "This does not prove the lithium problem is solved. It upgrades the evidence state from 'positive but ESS-limited' to 'CMB-gate audit-qualified'. Matched next-to-leading-order plasma physics, independent reproduction and an explicit particle/action-level model remain open.",
        ],
    )
    add_tv_index(
        doc,
        [
            ("Document update", "Regenerated 11 September 2026 from local post-process and combined qualification audit outputs."),
            ("Authors", "Clive Stewart Boyd; David Ng."),
            ("Computational support", "ChatGPT Codex sol 5.5 referenced as the project computational-assistance environment."),
            ("Primary new gate", "Full Planck PR4/NPIPE native nuisance audit: PASS."),
        ],
    )

    doc.add_heading("Current Claim Boundary", level=1)
    add_para(
        doc,
        "The current result should be described as an audit-qualified candidate for a selective electromagnetic pathway through the cosmological lithium problem. "
        "The analysis does not yet constitute a discovery claim, a completed microscopic model, or a replacement for standard BBN. "
        "It does show that the registered August candidate is not excluded by the native Planck PR4/NPIPE CMB nuisance posterior under the applied gates.",
    )
    add_table(
        doc,
        ["Evidence layer", "Status on 11 September 2026", "Interpretation"],
        [
            ["Standard BBN reproduction", "Pass", "The ordinary D/H and helium success, plus lithium excess, was reproduced."],
            ["July global-modification controls", "Fail for solution", "Useful negative controls; no D+He+Li pass."],
            ["Selective mass-7 threshold window", "Pass as sanity gate", "Photons near the Be-7 threshold can affect mass-7 more selectively than late neutrons."],
            ["ACROPOLIS/LINX/FIRAS candidate", "Pass at frozen candidate", "Abundance and distortion gates are satisfied in the recorded candidate cell."],
            ["Full Planck nuisance CMB gate", "Pass after ESS top-up", "Combined chain passes split, ESS, registered density and nuisance stability gates."],
            ["Matched NLO plasma calculation", "Open", "Needed before strong physics promotion."],
            ["Independent reproduction", "Open", "Needed before external paper-strength claim."],
        ],
        "Table 2. Evidence ladder status",
    )

    add_chronology(doc)

    doc.add_heading("Scientific Background", level=1)
    add_para(
        doc,
        "BBN predicts light-element abundances from the baryon density, the radiation content, neutron-proton conversion, nuclear reaction rates and the expansion history. "
        "The same framework that gives successful deuterium and helium predictions leaves a persistent Li-7 excess when compared with warm metal-poor halo-star plateaux. "
        "Most primordial Li-7 is produced first as Be-7, so a selective mechanism can target mass-7 after ordinary nucleosynthesis without needing to disrupt all light nuclei.",
    )
    doc.add_heading("Core Equations", level=2)
    add_table(
        doc,
        ["Label", "Expression", "Use"],
        [
            ["E1", "eta10 = 273.78 omega_b", "CMB-to-BBN baryon-density conversion for T_CMB = 2.7255 K."],
            ["E2", "dY_i/dt = sum_j N_ij Gamma_j(Y,T)", "Network abundance evolution in schematic form."],
            ["E3", "Li discrepancy factor = (Li-7/H)_BBN / (Li-7/H)_plateau", "Measures the anomaly scale."],
            ["E4", "E_gamma approximately m_phi / 2", "Two-photon decay energy for the mediator candidate."],
            ["E5", "Delta chi2 = (x - mu)^T C^-1 (x - mu)", "Registered-cell density test in the CMB posterior."],
            ["E6", "C = L L^T; C^-1 = L^-T L^-1", "SPD covariance cache used for the Planck PR4/NPIPE likelihood."],
        ],
        "Table 3. Formula cross-reference",
    )
    add_tv_index(
        doc,
        [
            ("Physics boundary", "Equations are diagnostic and model-selection support. A full Lagrangian/action is still required for theory-level promotion."),
            ("Symmetry boundary", "The BBN-Li mediator model is not currently claimed as a Noether-derived symmetry result."),
        ],
    )

    doc.add_heading("August Candidate Cell", level=1)
    add_para(
        doc,
        "The August candidate moved away from broad changes to the expansion law and toward a selective post-BBN transfer. "
        "Its working idea is simple: an electromagnetic injection can be tuned to interact with mass-7, principally Be-7, more efficiently than with deuterium, provided the spectrum and timing remain narrow enough.",
    )
    add_table(doc, ["Quantity", "Value"], candidate_rows(), "Table 4. ACROPOLIS/LINX candidate and CMB comparison")
    add_tv_index(
        doc,
        [
            ("Candidate source", "Post-process report and local resummed ACROPOLIS/LINX candidate values."),
            ("Abundance anchor", "Li-7/H = 1.488554e-10 in candidate units table, with D/H and Y_p retained inside gates."),
            ("Spectral-distortion anchor", "mu_distortion = 1.476434e-09, below the adopted FIRAS-style gate."),
        ],
    )

    doc.add_heading("Full Planck PR4/NPIPE CMB Nuisance Gate", level=1)
    add_para(
        doc,
        "The full CMB gate replaced the earlier compressed Gaussian approximation with a native spectrum-level Planck likelihood. "
        "It samples eta10 through omega_b, samples N_eff and Y_p independently, and includes the Planck high-l CamSpec TTTEEE likelihood, low-l TT, low-l EE and PR4 lensing, plus foreground and calibration nuisance parameters.",
    )
    add_table(
        doc,
        ["Component", "Role"],
        [
            ["Planck PR4/NPIPE CamSpec TTTEEE", "High-multipole CMB spectra and covariance."],
            ["Planck 2018 low-l TT", "Large-scale temperature likelihood."],
            ["Planck 2018 low-l EE", "Large-scale polarisation likelihood."],
            ["Planck PR4 lensing", "Lensing reconstruction constraint."],
            ["Foreground/calibration nuisance block", "Native amplitudes, tilts and calibration parameters left free under their likelihood priors."],
        ],
        "Table 5. Planck likelihood stack",
    )
    add_table(doc, ["Parameter", "Mean", "Std", "Median", "68% interval", "95% interval"], posterior_rows(), "Table 6. Main full CMB nuisance posterior summary")
    add_image_if_exists(
        doc,
        POST_ROOT / "planck_pr4_full_nuisance_key_posteriors.png",
        "Figure 1. Key posterior summaries from the Planck PR4/NPIPE nuisance run.",
    )
    add_tv_index(
        doc,
        [
            ("Protocol", str(OUTPUT_ROOT / "full_cmb_nuisance_protocol_2026-08-21.md")),
            ("Post-process", str(POST_ROOT / "planck_pr4_full_nuisance_postprocess_report.md")),
            ("Native-likelihood caveat", "This is not a Planck map-level reanalysis and does not include ACT or SPT."),
        ],
    )

    doc.add_heading("September Qualification Audit", level=1)
    add_para(
        doc,
        "The September audit combines the completed main chain with the ESS top-up chain. "
        "It uses a stricter gate set than the initial post-process: split-chain posterior summaries, autocorrelation-aware ESS, registered-cell posterior density, nuisance-parameter stability, all-row versus last-half comparison and explicit pass/fail scoring.",
    )
    add_table(doc, ["Gate", "Threshold", "Observed", "Pass"], gate_rows(), "Table 7. Preregistered qualification gates")
    add_table(doc, ["Parameter", "Split R-1", "Autocorr ESS", "tau_int", "Half-shift sigma", "Status"], metric_rows(), "Table 8. Selected audit metrics")
    add_table(doc, ["Quantity", "Value"], registered_density_rows(), "Table 9. Registered lithium cell in the combined CMB posterior")
    add_image_if_exists(
        doc,
        QUAL_ROOT / "planck_pr4_full_nuisance_audit_gates.png",
        "Figure 2. Qualification gate readout for the combined main plus ESS top-up chain.",
    )
    add_image_if_exists(
        doc,
        QUAL_ROOT / "planck_pr4_full_nuisance_split_posterior_overlay.png",
        "Figure 3. Split-chain posterior overlay used to audit stability.",
    )
    add_tv_index(
        doc,
        [
            ("Combined chain", str(OUTPUT_ROOT / "chains_combined" / "planck_pr4_full_nuisance_main_plus_ess_topup.1.txt")),
            ("Audit report", str(QUAL_ROOT / "planck_pr4_full_nuisance_qualification_audit.md")),
            ("Headline", "overall_pass=True; min autocorrelation ESS=1449; max split R-1=0.002143; registered Delta chi2=1.257."),
            ("Caveat", "Cobaya log R-1 line is inherited from the main run log; independent split R-1 and ESS are recomputed from the combined chain."),
        ],
    )

    doc.add_heading("Interpretation", level=1)
    add_para(
        doc,
        "The main interpretation has changed since 28 August 2026. At that stage the Planck nuisance gate was still pending. "
        "After the main chain completed on 2 September 2026 it was positive but not fully qualified because the stricter autocorrelation-aware ESS gate failed for some nuisance directions. "
        "The ESS top-up, combined with the main chain and audited on 11 September 2026, removes that immediate technical objection.",
    )
    add_para(
        doc,
        "The data therefore support the statement that the registered lithium candidate lies in an acceptable part of the native Planck PR4/NPIPE posterior for eta10, N_eff and Y_p. "
        "This is narrower than claiming the mediator exists. It is nevertheless materially stronger than the August status because it passes the spectrum-level CMB nuisance treatment rather than relying on a compressed CMB approximation.",
    )

    doc.add_heading("Future Research and Test Options", level=1)
    add_numbered(
        doc,
        [
            "Matched next-to-leading-order thermal calculation for the same mediator point, including pole-cut consistency and thermal self-energy treatment.",
            "Independent reproduction of the ACROPOLIS/LINX abundance chain and the Planck nuisance posterior audit on a separate machine or environment.",
            "Sensitivity analysis over nuclear-rate choices, stellar plateau assumptions, D/H observational anchors and helium determinations.",
            "Explicit action-level particle model: specify the field content, Lagrangian, production mechanism, decay channel, symmetries and conservation laws.",
            "External CMB and spectral tests: ACT/SPT consistency, improved FIRAS-style spectral-distortion treatment and future PIXIE-like sensitivity projection.",
            "Publication-grade reproducibility package: frozen manifests, file hashes, environment locks, CSV tables and figure regeneration commands.",
        ],
    )

    add_common_appendices(doc)
    add_footer(doc)
    doc.save(REPORT_PATH)


def build_formal_paper() -> None:
    doc = setup_doc("Formal scientific paper BBN lithium August findings 2026-09-11")
    add_title_block(
        doc,
        "A Selective Electromagnetic-Transfer Candidate for the Cosmological Lithium Problem",
        "August 2026 findings with September 2026 full Planck PR4/NPIPE nuisance qualification",
        "Formal paper draft. CMB gate qualified; mechanism promotion remains conditional.",
    )

    doc.add_heading("Lay Person Description", level=1)
    add_para(
        doc,
        "The early Universe made mostly hydrogen and helium, plus tiny amounts of deuterium, helium-3 and lithium. "
        "Modern calculations predict deuterium and helium very well, but they predict too much lithium compared with old metal-poor stars. "
        "This long-standing mismatch is the cosmological lithium problem.",
    )
    add_para(
        doc,
        "The July tests tried broad ways to change the early-Universe calculation. Those did not work: when lithium moved in the right direction, something else, usually deuterium or helium, moved out of agreement. "
        "The August work then tested a more selective idea: a late, narrow electromagnetic process could reduce the mass-7 reservoir after normal nucleosynthesis without spoiling the successful deuterium and helium predictions.",
    )
    add_para(
        doc,
        "The September update asks whether this selected candidate is still acceptable when compared with the detailed Planck CMB data, including calibration and foreground nuisance parameters. "
        "The answer is yes under the current audit gates. The result is not yet proof of a new particle or interaction, but it is strong enough to justify a formal follow-up programme.",
    )

    doc.add_heading("Connection to Existing Papers", level=1)
    add_para(
        doc,
        "The work sits inside a long literature on the lithium anomaly, precision BBN, stellar depletion, electromagnetic cascades, CMB constraints and spectral distortions. "
        "The references listed below are used as the current interpretive frame. They define the observational problem, the standard BBN baseline, the cascade machinery and the CMB/spectral gates.",
    )
    add_table(
        doc,
        ["Topic", "Representative references", "Role in this paper"],
        [
            ["Problem definition", "Fields 2011; Cyburt, Fields and Olive 2008; Spite and Spite 1982; Sbordone et al. 2010; Ryan 2023.", "Defines the Li-7 excess relative to metal-poor stellar plateaux."],
            ["Precision BBN", "Cyburt et al. 2016; Pitrou et al. 2018; Fields et al. 2020; LINX 2024.", "Defines standard abundance calculations and CMB-compatible baryon density handling."],
            ["Observational anchors", "Cooke et al.; helium determinations; Planck Collaboration 2020.", "Defines D/H, Y_p and CMB posterior constraints."],
            ["Cascade and late-decay physics", "Depta et al. ACROPOLIS 2021 and related electromagnetic-decay constraints.", "Tests whether a late electromagnetic process damages other light elements."],
            ["CMB and spectral gates", "Planck PR4/NPIPE/CamSpec; Fixsen FIRAS; Khatri and Sunyaev.", "Tests CMB consistency and spectral-distortion safety."],
        ],
        "Table 1. Literature cross-reference",
    )

    doc.add_heading("Abstract", level=1)
    add_para(
        doc,
        "A two-stage investigation of the cosmological lithium problem is reported. The July 2026 phase tested standard and modified BBN controls, including expansion, clock-proxy, reaction-rate and SU2-style variations, and found no joint D/H, helium and Li-7 solution. "
        "The August 2026 phase instead tested a selective post-BBN electromagnetic-transfer candidate. A mediator cell with m_phi = 4.44 MeV, E_gamma = 2.22 MeV, tau_phi = 1.584893e5 s and n_phi/n_gamma = 1.211528e-6 gives Li-7/H = 1.488554e-10 while retaining the recorded D/H, Y_p, He3/D, entropy and spectral-distortion gates. "
        "A native Planck PR4/NPIPE CMB nuisance run subsequently sampled omega_b, N_eff and Y_p independently with the CamSpec TTTEEE, low-l TT, low-l EE and PR4 lensing likelihoods. "
        "The combined main plus ESS top-up audit passes all preregistered qualification gates, including min autocorrelation-aware ESS = 1449, independent split R-1 max = 0.002143 and registered-cell Delta chi2 = 1.257. "
        "The candidate is therefore CMB-gate qualified but remains conditional on matched NLO thermal physics and independent reproduction.",
    )

    doc.add_heading("1. Introduction", level=1)
    add_para(
        doc,
        "The cosmological lithium problem is unusual because standard BBN is not globally broken. Deuterium and helium agree broadly with a CMB-inferred baryon density, yet Li-7/H remains high relative to the stellar plateau. "
        "Any proposed mechanism must therefore be selective. A broad change to expansion or baryon density is heavily constrained by the same observables that standard BBN already explains.",
    )
    add_para(
        doc,
        "This paper records the transition from negative July controls to a positive August candidate and then to a September CMB-nuisance qualification. "
        "The chronology matters because it prevents the August candidate being mistaken for an arbitrary retuning of the failed July pathways.",
    )

    doc.add_heading("2. Methods", level=1)
    add_para(
        doc,
        "The computational programme used LINX for BBN network and abundance validation, ACROPOLIS for electromagnetic cascade checks, Cobaya and CAMB for native Planck CMB likelihood sampling, and project-specific Python audit layers for qualification gates and reproducibility. "
        "The analysis separates three claims: reproducing the standard lithium problem, finding a parameterised selective transfer that repairs lithium without damaging other gates, and promoting that parameterisation to a microscopic theory.",
    )
    add_table(
        doc,
        ["Method block", "Implementation", "Decision role"],
        [
            ["Standard and modified BBN controls", "July LINX/FR/SU2 scripts.", "Negative control and failure-boundary definition."],
            ["Selective transfer/cascade checks", "ACROPOLIS plus LINX candidate pipeline.", "Candidate discovery and abundance safety screen."],
            ["Finite-temperature production", "Momentum-resolved and resummed production calculations.", "Tests whether the candidate abundance is physically plausible."],
            ["Native CMB nuisance likelihood", "Cobaya/CAMB with Planck PR4/NPIPE CamSpec, low-l TT, low-l EE and PR4 lensing.", "Replaces compressed CMB gate."],
            ["Qualification audit", "Split-chain, autocorrelation ESS, nuisance stability and registered-cell density scripts.", "Determines whether the CMB gate is technically qualified."],
        ],
        "Table 2. Method blocks",
    )

    doc.add_heading("3. Results", level=1)
    doc.add_heading("3.1 July Negative Controls", level=2)
    add_para(
        doc,
        "The July programme reproduced the anomaly and tested whether broad modifications could resolve it. The result was negative: no tested global modification simultaneously preserved deuterium, helium and lithium. "
        "This result is scientifically useful because it sharply constrains what a successful mechanism can look like.",
    )
    add_tv_index(
        doc,
        [
            ("Stage", "July 16-18, 2026."),
            ("Result", "No D/H + helium + Li-7 pass from standard, FR/clock-proxy, reaction-rate or SU2-style control scans."),
            ("Consequence", "Future solution must be selective rather than a broad background retuning."),
        ],
    )

    doc.add_heading("3.2 August Selective Candidate", level=2)
    add_para(
        doc,
        "The August candidate uses a narrow electromagnetic window after ordinary BBN. Its central physical premise is that mass-7, mostly Be-7 before electron capture, can be reduced without crossing into broad deuterium destruction. "
        "The candidate therefore targets the unusual topology of the mass-7 pathway rather than changing the whole thermal history.",
    )
    add_table(doc, ["Quantity", "Value"], candidate_rows(), "Table 3. Registered August candidate")
    add_tv_index(
        doc,
        [
            ("Stage", "August 20-21, 2026."),
            ("Positive result", "Candidate passes the recorded abundance/cascade/FIRAS gates at the frozen cell."),
            ("Boundary", "A complete action-level particle model is still not supplied."),
        ],
    )

    doc.add_heading("3.3 Native Planck CMB Nuisance Qualification", level=2)
    add_para(
        doc,
        "The native Planck PR4/NPIPE test is the strongest new constraint in this update. It lets omega_b, N_eff and Y_p vary independently, and marginalises over foreground and calibration nuisance directions rather than using only a compressed two-variable Planck prior.",
    )
    add_table(doc, ["Gate", "Threshold", "Observed", "Pass"], gate_rows(), "Table 4. Qualification gates")
    add_table(doc, ["Quantity", "Value"], registered_density_rows(), "Table 5. Registered-cell posterior density")
    add_table(doc, ["Parameter", "Mean", "Std", "Median", "68% interval", "95% interval"], posterior_rows(), "Table 6. CMB posterior summary")
    add_image_if_exists(
        doc,
        QUAL_ROOT / "planck_pr4_full_nuisance_audit_gates.png",
        "Figure 1. Combined main plus ESS top-up qualification gates.",
    )
    add_tv_index(
        doc,
        [
            ("Stage", "September 11, 2026."),
            ("Combined chain size", "120073 compressed rows; raw weight sum 393272."),
            ("Gate outcome", "PASS on all preregistered gates."),
            ("Most conservative ESS", "1449 for amp_143x217 among audited gate parameters."),
        ],
    )

    doc.add_heading("4. Discussion", level=1)
    add_para(
        doc,
        "The September audit changes the evidence state but not the claim boundary. It demonstrates that the registered lithium candidate is compatible with the detailed Planck CMB nuisance posterior under the applied gates. "
        "It does not show that the mediator is realised in nature, nor that the calculated transfer is stable under every possible plasma correction.",
    )
    add_para(
        doc,
        "The candidate is promising because it respects the selectivity demanded by the lithium problem. It is risky because narrow electromagnetic windows and late-decay mechanisms are sensitive to cascade modelling, abundance priors, stellar interpretation and spectral-distortion assumptions. "
        "Those risks are manageable but must be handled before publication as a strong physical claim.",
    )

    doc.add_heading("5. Limitations", level=1)
    add_bullets(
        doc,
        [
            "The analysis is not yet an independent external reproduction.",
            "The mediator model is still phenomenological; a fully specified action, symmetry content and particle-production history remain required.",
            "The Planck run is spectrum-level but not map-level, and does not include ACT or SPT spectra.",
            "The abundance and CMB gates are joined by a registered cell rather than by a single global sampler over all microscopic nuisance directions.",
            "Systematic uncertainty in the stellar lithium plateau remains a separate astrophysical issue.",
        ],
    )

    doc.add_heading("6. Conclusions", level=1)
    add_para(
        doc,
        "The July investigations failed to resolve the lithium problem using broad BBN or cosmological modifications. The August investigations found a selective electromagnetic-transfer candidate. "
        "The September full Planck PR4/NPIPE nuisance audit now qualifies the CMB gate for that candidate with all preregistered audit tests passed. "
        "The correct conclusion is that the candidate merits formal follow-up, not that the lithium problem is closed.",
    )

    doc.add_heading("7. Future Work", level=1)
    add_numbered(
        doc,
        [
            "Run matched next-to-leading-order thermal corrections at and around the registered candidate cell.",
            "Reproduce the ACROPOLIS/LINX candidate with independent code and fixed public manifests.",
            "Extend CMB checks to ACT/SPT where likelihood access and parameter consistency allow.",
            "Publish a frozen reproducibility bundle with checksums for scripts, input tables, plots and DOCX/PDF outputs.",
            "Promote the model only after an action-level declaration specifies fields, interactions, gauge content, conservation laws and priors.",
        ],
    )

    add_common_appendices(doc)
    add_footer(doc)
    doc.save(PAPER_PATH)


def copy_to_github() -> None:
    if not GITHUB_PAPERS.exists():
        return
    for path in [REPORT_PATH, PAPER_PATH]:
        if path.exists():
            shutil.copy2(path, GITHUB_PAPERS / path.name)
    if GITHUB_CODE.exists():
        shutil.copy2(Path(__file__).resolve(), GITHUB_CODE / Path(__file__).name)


def main() -> None:
    build_update_report()
    build_formal_paper()
    copy_to_github()
    print(f"Saved update report: {REPORT_PATH}")
    print(f"Saved formal paper: {PAPER_PATH}")
    if GITHUB_PAPERS.exists():
        print(f"Copied DOCX files to: {GITHUB_PAPERS}")
    if GITHUB_CODE.exists():
        print(f"Copied builder to: {GITHUB_CODE / Path(__file__).name}")


if __name__ == "__main__":
    main()
