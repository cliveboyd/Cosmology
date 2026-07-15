# Research Paper Format Profile

Reference document: `SU2_Quaia_Chiral_Gauge_Sector_Technical_Paper_2026-07-15.docx`

Use this as the default format for future research-paper regenerations unless the requested paper has a stronger existing template.

## Structure

- Start with a compact title page and metadata table.
- Put a clear `Claim Boundary` block before the abstract.
- Use the abstract as a short technical summary, not as the first place where caveats appear.
- Keep the hierarchy simple: `Heading 1`, `Heading 2`, `Heading 3`, `Normal`, and tables.
- Prefer one main topic per page or half-page region where practical.
- Put interpretation text immediately after the table it explains.
- End with reproducibility/archive artifacts, source links, and acceptance gates.

## Tone

- Shorter paragraphs than the generated draft.
- Direct claim language: distinguish target-follow-up, control, caveat, and promotion gate.
- Avoid discovery phrasing unless the global and systematics gates actually pass.
- Keep equations and model definitions concise; do not let formulas dominate the page.

## Layout

- Restrained formal report style with blue section headings and quiet gray tables.
- Tables are the evidence carrier; prose explains what the table means.
- Use shaded callout/interpretation rows sparingly, mainly for claim boundaries and bottom-line readouts.
- Avoid orphaned table rows at page starts.
- Avoid blank spacer pages.
- Keep formula boxes left-aligned and compact.

## Reproducibility

- Include local artifact paths, program/script inventory, datasets, and manifest references.
- Keep citation-facing PDFs in `github_export/docs/papers/`.
- Keep editable DOCX sources in the Codex outputs folder unless explicitly asked to track them.
- Record SHA-256 hashes and page counts in `github_export/manifests/`.
