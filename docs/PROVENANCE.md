**Source identity and changes**

Prepared on 18 September 2026 from the recovered project artifacts.

| Component | Identity |
|---|---|
| Bundled `vendor/eT-development.zip` | SHA-256 `45724c57cf4f3c28c889e0749f650d9aaef49e3c99690361525073b458dbe92d` |
| Original speed-v2 workflow | Archive SHA-256 `499575ff3af0a95f4ecefe8c5c69396ddb3b2b05ed32498e08592c7c02f6da43` |
| Original Eckart recovery archive | SHA-256 `d84e137d0530f483412c534868c3c1eb5c1322d55341d0aac6cc24db22744b4b` |
| i-PI | Python distribution `ipi==3.2.0`; [official release](https://github.com/i-pi/i-pi/releases/tag/v3.2.0) |
| libcint | [Upstream tag v6.1.2](https://github.com/sunqm/libcint/releases/tag/v6.1.2); actual checkout commit recorded by installer |
| LibXC, optional | [Upstream tag 5.1.3](https://gitlab.com/libxc/libxc/-/tree/5.1.3); actual checkout commit recorded by installer |
| Python and build tooling | Direct versions in `requirements.txt`; resolved environment recorded with `pip freeze` at installation |

The supplied eT archive is the recovered development snapshot, identified by its bytes, not asserted to be the latest upstream commit. The source archive remains unmodified in this distribution. Installation extracts it and applies the recovered fixed-coherent-state gradient patch. No new electronic-gradient or DFT implementation was added. eT's author notices and GPLv3 license are retained.

The workflow starts from speed-v2's SCF restart, fallback, rolling scratch and parser code. This revision adds:

- Direct execution of the native eT binary with explicit basis and OpenMP environment settings.
- Cold fallback inputs that omit the restart directive; validated restart files are promoted only after output parsing.
- Correct mixed-unit treatment of standalone XYZ files: nuclei in angstrom, final L.x in photon atomic units.
- Normal handling of i-PI's idle socket closure, while truncated protocol messages remain errors.
- A supervisor that derives symbols and the photon index from embedded bead arrays, validates photon mass and unused DOFs, runs force checks and supervises its own process groups.
- Global Eckart preparation for any nonlinear nuclear group, including 21 waters, with mass-weighted momentum projection and an explicit fixed reference preserved from existing Eckart checkpoints.
- Native/launcher comparisons, electronic force validation, thread benchmarking and NVE drift summaries.
- A guard against treating ordinary DFT inputs as a validated CBO force provider.

The latest recovered `i_pi_eckart_compat.py` is included unchanged. It supplies in-memory fixes for stale fixed-DOF parameter names and the MTS force accessor in i-PI's constrained path. This bundle pins and tests i-PI 3.2.0; other versions are not certified by statements in older recovery documentation.

The original monomer-specific Eckart preparation script and its recovery README are retained under `workflow/legacy_eckart/`. Their historical paths and validation statements describe the recovery artifact, not the new installer. Use the top-level README and `workflow/prepare_run.py` for the current workflow.

`workflow/examples/water21.xyz` is the user's recovered 21-water configuration with a final L particle. Its inclusion does not certify equilibration, stability or suitability for a specific thermodynamic state. The generated XML examples contain reproducible random thermal momenta and are examples rather than equilibrated production restarts.

The installer uses eT's documented [OpenMP and BLAS integer settings](https://etprogram.org/setup.html) and the matching source's CMake interface. libcint uses the ordering and screening compile options required by eT's [integral-library instructions](https://etprogram.org/integral_libraries.html). The CPU build uses optimization without introducing a new fast-math approximation.

The [official DFT documentation](https://etprogram.org/dft_calculation.html) documents ordinary restricted single-point DFT. No recovered source branch or artifact established a complete CBO-DFT force implementation. Linking LibXC alone therefore does not enable the requested future cavity-DFT MD route.
