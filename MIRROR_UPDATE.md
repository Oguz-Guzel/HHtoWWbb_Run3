# How to Update the Mirrored Repository

This repository on GitHub (`https://github.com/Oguz-Guzel/HHtoWWbb_Run3`) is a mirror of the original GitLab repository (`ssh://git@gitlab.cern.ch:7999/aguzel/HHtoWWbb_Run3.git`).

## Prerequisites

Before updating the mirror, ensure you have:
- Access to both the GitLab (original) and GitHub (mirror) repositories
- Git installed and configured on your local machine
- SSH keys set up for both GitLab CERN and GitHub

## Setup Mirror Repository (One-time Setup)

If you're setting up the mirror for the first time:

```sh
# Clone the repository with both remotes
git clone https://github.com/Oguz-Guzel/HHtoWWbb_Run3.git
cd HHtoWWbb_Run3

# Add the GitLab repository as upstream remote
git remote add upstream ssh://git@gitlab.cern.ch:7999/aguzel/HHtoWWbb_Run3.git

# Verify remotes
git remote -v
```

You should see:
```
origin    https://github.com/Oguz-Guzel/HHtoWWbb_Run3 (fetch)
origin    https://github.com/Oguz-Guzel/HHtoWWbb_Run3 (push)
upstream  ssh://git@gitlab.cern.ch:7999/aguzel/HHtoWWbb_Run3.git (fetch)
upstream  ssh://git@gitlab.cern.ch:7999/aguzel/HHtoWWbb_Run3.git (push)
```

## Updating the Mirror

To sync the latest changes from GitLab to GitHub:

### Method 1: Update Specific Branch (Recommended)

```sh
# Make sure you're on the main branch (or the branch you want to update)
git checkout main

# Fetch latest changes from GitLab
git fetch upstream

# Merge upstream changes into your local branch
git merge upstream/main

# Push the updated branch to GitHub
git push origin main
```

### Method 2: Update All Branches and Tags

```sh
# Fetch all branches and tags from GitLab
git fetch upstream --tags

# For each branch you want to mirror, checkout and update
git checkout main
git merge upstream/main
git push origin main

# Push all tags
git push origin --tags
```

### Method 3: Complete Mirror Sync (All Branches)

For a complete sync of all branches:

```sh
# Fetch everything from upstream
git fetch upstream

# Get list of all remote branches from upstream
git branch -r | grep upstream/ | grep -v '\->' | sed 's/upstream\///' | while read branch; do
    # Checkout or create local tracking branch
    git checkout -B "$branch" "upstream/$branch"
    # Push to GitHub
    git push -f origin "$branch"
done

# Push all tags
git push origin --tags --force

# Return to main branch
git checkout main
```

## Updating Git Submodules

This repository contains the `TransformerNN` submodule. To initialize and update it:

```sh
# Initialize submodules (first time only)
git submodule init

# Update submodules to the latest commit referenced in the repository
git submodule update

# Or do both in one command
git submodule update --init --recursive
```

To update the submodule to its latest version:

```sh
# Navigate to the submodule directory
cd TransformerNN

# Fetch and checkout the latest changes
git fetch origin
git checkout main  # or master, depending on the default branch
git pull origin main

# Go back to the main repository
cd ..

# Stage the submodule update
git add TransformerNN

# Commit the update
git commit -m "Update TransformerNN submodule to latest version"

# Push to both remotes
git push origin main
git push upstream main
```

## Automated Mirroring (Optional)

For automated mirroring, you can use GitHub Actions. Create `.github/workflows/mirror.yml`:

```yaml
name: Mirror GitLab to GitHub

on:
  schedule:
    # Run every 6 hours
    - cron: '0 */6 * * *'
  workflow_dispatch:  # Allow manual triggering

jobs:
  mirror:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout
        uses: actions/checkout@v3
        with:
          fetch-depth: 0
          
      - name: Configure Git
        run: |
          git config --global user.name "GitHub Action"
          git config --global user.email "action@github.com"
          
      - name: Add upstream remote
        run: |
          git remote add upstream ssh://git@gitlab.cern.ch:7999/aguzel/HHtoWWbb_Run3.git || true
          
      - name: Fetch from upstream
        run: git fetch upstream --tags
        
      - name: Merge and push
        run: |
          git checkout main
          git merge upstream/main
          git push origin main
          git push origin --tags
```

## Troubleshooting

### Authentication Issues

If you encounter authentication problems:

```sh
# For GitLab CERN, ensure your SSH key is properly configured
ssh -T git@gitlab.cern.ch

# For GitHub
ssh -T git@github.com
```

### Merge Conflicts

If there are conflicts when merging from upstream:

```sh
# Abort the merge
git merge --abort

# Alternatively, resolve conflicts manually
git status  # See conflicting files
# Edit the files to resolve conflicts
git add <resolved-files>
git commit
git push origin main
```

### Force Update (Use with Caution)

If the GitHub mirror has diverged and you want to make it identical to GitLab:

```sh
git fetch upstream
git checkout main
git reset --hard upstream/main
git push origin main --force
```

⚠️ **Warning**: Force pushing will overwrite the GitHub repository's history. Only use this if you're certain the GitLab repository is the source of truth.

## Best Practices

1. **Regular Updates**: Update the mirror regularly to keep it in sync
2. **Branch Protection**: Consider enabling branch protection on GitHub to prevent accidental direct commits
3. **Communication**: Inform collaborators which repository (GitLab or GitHub) is the primary source
4. **Submodule Tracking**: Keep submodules updated separately and commit those updates
5. **Tags**: Don't forget to push tags when mirroring releases

## Verification

After updating, verify the sync:

```sh
# Check that both remotes have the same commits
git log origin/main --oneline -5
git log upstream/main --oneline -5

# Check that they point to the same commit
git rev-parse origin/main
git rev-parse upstream/main
```

Both commands should return the same commit hash.

## Quick Reference

```sh
# Quick update workflow
git checkout main
git fetch upstream
git merge upstream/main
git push origin main
git push origin --tags
git submodule update --init --recursive
```
