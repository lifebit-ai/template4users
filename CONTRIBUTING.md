# Contributing Guidelines

## Bug-fix/Feature process

* First of all, please create an issue if it does not exist already. This will help create the visibility among the team to discuss further.
* Always use the `dev` branch to checkout from for a feature/bug-fix and then create PR against it.

## Git branching process

* `main` branch is only reserved for stable versions, so it should be only used for release purpose.
* Most bug fixes and features should be merged only onto `dev`. This will be implemented by opening PR reviews against this branch.

## Nextflow-specific process

Some specifications to ensure structure standards across Nextflow pipelines in the organisation.

* Tool/software dependencies of a pipeline are handled by combination of conda and docker (`environment.yml` + `Dockerfile`) in `containers/` directory
* Try to make the tool versions very specific in `environment.yml` file. Ex - `samtools==1.12`
* Try to add test profiles for each of the features making into the pipeline and also make sure they are part of github-actions Continues-Integration (CI) test (`.github/workflows/ci.yml`)
* Test profiles will need small datasets, as they quickly need to be tested in CI. If the dataset is less than few KB keep them under `testdata` in same repo as pipeline or move them to appropriate places in `s3://lifebit-featured-datasets/pipelines/`

## Release process

Checklist

* Version bump in:
  * In nextflow.config
  * In VERSION file
  * Change the tag for containers used (should be same with the VERSION file)
* Add release note to the changelog.md

NOTE: The docker containers with VERSION will be automatically built and pushed to quay.io
