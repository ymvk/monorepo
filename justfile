# justfile -- reproducible git + release ops for this monorepo.
# Requires: just, git, gh (GitHub CLI).

# Edit these if the remote owner, SSH alias, or repo name changes.
owner  := "ymvk"
remote := "github-personal"
repo   := `basename "$PWD"`

# Default: list available recipes.
default:
    @just --list

# Show who we are and where we push. Run before anything risky.
whoami:
    @echo "user:    $(git config user.name) <$(git config user.email)>"
    @echo "remote:  $(git remote get-url origin 2>/dev/null || echo '(none)')"

# First-time setup. Assumes github.com/{{owner}}/{{repo}} already exists on
# GitHub (create it in the UI first). Uses the {{remote}} SSH alias so both
# push auth and commit identity resolve to the second account -- identity
# comes from the conditional include in ~/.gitconfig, keyed on this alias.
init:
    @test ! -d .git || (echo "already a git repo" && exit 1)
    git init
    git branch -M main
    git remote add origin git@{{remote}}:{{owner}}/{{repo}}.git
    git add .
    git commit -m "first commit"
    git push -u origin main
    @just whoami

# Stage everything, commit (skipped if nothing staged), push.
push msg="update":
    @just whoami
    git add -A
    @if git diff --cached --quiet; then \
        echo "nothing staged"; \
    else \
        git commit -m "{{msg}}"; \
    fi
    git push

# Rewrite the last commit under the currently configured identity. Use once
# if a prior commit was attributed to the wrong account, then `push-force`.
amend-author:
    git commit --amend --reset-author --no-edit
    @just whoami

# Force-push with lease protection (safer than --force).
push-force:
    git push --force-with-lease

# Trigger the validate workflow manually. Leave `package` empty to validate
# all packages, or pass a folder name to scope to one.
validate package="":
    @if [ -z "{{package}}" ]; then \
        gh workflow run validate.yml; \
    else \
        gh workflow run validate.yml -f package="{{package}}"; \
    fi
    @sleep 3
    gh run watch --exit-status

# Cut a release tag in `<pkg>-v<version>` form and push it. A matching
# release workflow (release-pypi, release-ghcr) fires if the package has
# the relevant manifest or Dockerfile.
tag package version:
    @test -d "packages/{{package}}" || (echo "no such package: packages/{{package}}" && exit 1)
    git tag "{{package}}-v{{version}}"
    git push origin "{{package}}-v{{version}}"

# List existing tags for one package.
tags package:
    @git tag --list "{{package}}-v*"
