# Publish the repository and GitHub Pages

Suggested repository name: **cavity-qed-md**. Suggested description:

> eT–i-PI interface for cavity Born–Oppenheimer molecular dynamics and joint nuclear–photon path-integral sampling.

## Before a formal release

1. Compare the preserved source snapshot with the current working workstation directory, particularly any later eT patches and PLUMED files.
2. Select and record the interface license. Preserve all upstream licenses and author notices. Citation metadata for the accompanying manuscript is already recorded in `CITATION.cff`, the README and the website; add its permanent ChemRxiv identifier after posting.
3. Attach real-eT force-validation results and representative NVE/PIMD logs to the version they validate.

The files can be reviewed privately before choosing public visibility. This package does not include a made-up public URL or claim a tagged production release.

## Create and push

Create an empty repository under the intended GitHub user or organization. From the extracted project root:

```bash
git init -b main
git add .
git commit -m "Prepare cavity-QED MD source and documentation"
git remote add origin https://github.com/OWNER/cavity-qed-md.git
git push -u origin main
```

Replace `OWNER` with the intended account. Do not initialize the remote with a second README if pushing this tree into a new empty repository.

## Enable Pages

In the repository, open **Settings → Pages → Build and deployment**. Select **Deploy from a branch**, branch **main**, folder **/docs**, then save. The website uses plain HTML, inline CSS/JavaScript and `.nojekyll`; no Node build or external font service is required.

The website's source download is served from `docs/downloads/cavity-qed-md-source.tar.gz`. SHA-256 verification is available beside it. After Pages reports success, use GitHub's returned URL; do not infer that publication succeeded from the source push alone.

Official instructions: https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site

## Maintain the download

After posting the manuscript on ChemRxiv, add its verified DOI and public URL to `CITATION.cff` (under `preferred-citation`), `docs/cavity-qed-md.bib`, the README citation section and the website's `#citation` section. Update the unpublished status and BibTeX entry consistently with the ChemRxiv record, preserving the published title and author order.

After editing source or documentation, run:

```bash
python3 scripts/build_release.py
python3 scripts/check_release.py
```

Commit the refreshed source archive and manifests along with the edits. The script excludes run output, credentials, caches and Git metadata. The GitHub Actions check runs the existing interface tests and verifies the prepared files; it does not build eT or certify electronic forces.
