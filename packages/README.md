# packages

## Overview

A package is any folder directly under `packages/`. Language is detected from
its manifest file:

| File              | Language |
| ----------------- | -------- |
| `pyproject.toml`  | Python   |
| `Cargo.toml`      | Rust     |
| `build.zig`       | Zig      |
| `dune-project`    | OCaml    |
| `DESCRIPTION`     | R        |
| `CMakeLists.txt`  | C++      |
| `Makefile`        | C        |

Packages are discovered via CI on every push and runs the language-appropriate test
command. See `.github/workflows/ci.yml`.

`.github/workflows/validate.yml` exercises the full publish path without
pushing anything: `python -m build` + `twine check --strict` for Python,
`cargo publish --dry-run` for Rust, and `podman build` + `podman save` for
any package with a `Dockerfile`. It catches broken manifests, missing files
in the sdist, unrenderable READMEs, and Dockerfile regressions before a real
release tag.

`validate` is opt-in — it never runs on push or PR, only when you trigger
it manually. This keeps idle Actions usage at zero.

### Running validate

From the GitHub UI:

1. Go to the **Actions** tab.
2. In the left sidebar, click **validate**.
3. Click **Run workflow** (top-right of the run list).
4. *(Optional)* Fill in the **package** input with a single package name
   (e.g. `foo`) to validate just that one. Leave empty to validate all.
5. Click **Run workflow** to confirm.

From the command line (requires [`gh`](https://cli.github.com)):

```bash
# validate everything
gh workflow run validate.yml

# validate one package
gh workflow run validate.yml -f package=foo

# watch the run
gh run watch
```

### Using a package

```bash
# Python
pip install "foo @ git+https://github.com/<owner>/monorepo@<ref>#subdirectory=packages/foo"

# Rust (in Cargo.toml)
foo = { git = "https://github.com/<owner>/monorepo", rev = "<sha>" }

# R
remotes::install_github("<owner>/monorepo", subdir = "packages/foo", ref = "<ref>")

# OCaml
opam pin add foo https://github.com/<owner>/monorepo.git#<ref>
```

### Publishing a release

Tags follow the format `<package>-v<version>`, e.g. `foo-v1.0.0`. A matching tag triggers the relevant release workflow:

- `.github/workflows/release-pypi.yml` — publishes Python packages to PyPI via
  trusted publishing (configure the project on PyPI once, no tokens needed).
- `.github/workflows/release-ghcr.yml` — builds and pushes a container image for any package that has a `Dockerfile`. Consume with `podman pull ghcr.io/<owner>/<pkg>:<version>`.

Tags are optional. Install-from-git by commit SHA works without them.
