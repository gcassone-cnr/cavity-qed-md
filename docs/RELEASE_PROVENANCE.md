# Release preparation and source identity

Prepared on 30 October 2026 for the project **cavity-QED Molecular Dynamics**.

## Recovered inputs

| Original artifact | SHA-256 |
|---|---|
| `et_ipi_workstation_20260918.tar.gz` | `bcd378da3055ff99e27c4e2ddf706c3c2db39e5789217866fb17ffd853795cb3` |
| `h2o-cbohf-pimd.tar.gz` | `1cf47f79ce01c11cec3792441665163305fba826f1aa84834c78ccf8921f73ab` |
| Bundled `vendor/eT-development.zip` | `45724c57cf4f3c28c889e0749f650d9aaef49e3c99690361525073b458dbe92d` |

The workstation bundle supplies `install.sh`, `requirements.txt`, `workflow/`, `vendor/`, the original `docs/PROVENANCE.md`, and the original validation artifacts. Its README is retained verbatim as `ORIGINAL_INSTALLATION.md`; `docs/INSTALLATION.md` adapts its extraction instructions and documentation paths for this prepared source tree.

The PIMD bundle supplies the separate `examples/monomer-pimd/` tree. This retains the historical driver and launcher rather than silently merging them with the later workstation workflow.

The driver, scientific analysis code, original tests, templates, examples and vendored eT source archive retain their recovered bytes. New files add a landing page, documentation, publishing instructions, CI and packaging checks. The installer was subsequently corrected to run eT's own source generator before direct CMake configuration, following a reported missing `eri_tool_c_class.F90` error. The fixed-coherent gradient patch and its application script were then upgraded to v4 after a real-eT force check exposed missing fixed-displacement DSE derivatives; the same repair is included in both preserved driver distributions. See `GRADIENT_FIX_V4.md` for evidence and scope. The original workstation manifest is retained as `validation/ORIGINAL_MANIFEST.sha256`; the new root `SHA256SUMS` covers the distributable files.

These preserved snapshots may precede files currently used on the workstation. In particular, the current PLUMED integration is not present. An eventual tagged release should identify the exact workstation revision it reproduces.

See the original [source provenance](PROVENANCE.md) for pinned dependency identities and the history of the fixed-coherent gradient patch.
