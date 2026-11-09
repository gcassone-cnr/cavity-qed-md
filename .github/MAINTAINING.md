# Maintaining the documentation and source distribution

GitHub Pages serves `docs/` from the `main` branch. The site uses static HTML, CSS and JavaScript with `.nojekyll`.

After changing source or documentation, refresh and verify the distribution from the repository root:

```bash
python3 scripts/build_release.py
python3 scripts/check_release.py
```

Commit the refreshed archive and both checksum manifests with the changes. The source download is `docs/downloads/cavity-qed-md-source.tar.gz`. Verify the live site and download after GitHub Pages reports a successful deployment.

Keep citation metadata consistent across `CITATION.cff`, `docs/cavity-qed-md.bib`, the README and the website. Preserve upstream license texts and author notices. Record source and executable identities with validation results.

The GitHub Actions workflow checks the source package and interface. Use the compiled eT executable for the electronic force checks in `docs/VALIDATION.md`.
