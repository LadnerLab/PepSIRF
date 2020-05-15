# Contributing -- Steps for Contributing to PepSIRF and Releasing New Versions

# Introduction 
This document details the steps that are necessary to contribute to PepSIRF. 
This includes an outline of the git strategy that is used 
([GitFlow](https://nvie.com/posts/a-successful-git-branching-model/)), 
commands that are necessary, etc. 
If the steps in this document are followed closely, current and future members of 
LadnerLab (and the general public, if applicable) should be able to contribute to PepSIRF 
in a manner that allows many concurrent contributions, while maintaining 
reproducibility via [semantic versioning](https://semver.org/).

# GitFlow 
This project uses [GitFlow](https://nvie.com/posts/a-successful-git-branching-model/) 
as a strategy for effective use of Git. 
It is important to familiarize yourself with the reasoning behind GitFlow's 
strategy, and to be comfortable with the idea of the branching structure that 
is used before attempting to contribute to PepSIRF. 

# Semantic Versioning 
[Semantic Versioning](https://semver.org/) (semver) is a semantic strategy to 
handle project version numbers, and defines rules that dictate how the version number 
of a project should change as the project evolves. PepSIRF uses semver to handle the 
version information of any updates that are released. 

# Developing a new feature/bug fix
When working on a new feature or bug fix, it is important that a new child 
branch of ```develop``` is created. For naming of this new branch, we introduce the 
following notation:

- For branches that will implement a new feature, the ```feature-XXX``` format 
should be used, where ```XXX``` is a very brief description of the feature.
An example is ```feature-demux_input_verify```, which is one possible name for a 
branch that will implement a feature for the demux module that handles verification of 
input files.

- For branches that will implement a bug fix, the ```bug-XXX``` format should 
be used, where ```XXX``` is a very brief description of the bug. 
For example, ```bug-subjoin_mem_usage``` can be used for a bug that is 
related to the subjoin module's memory usage. 

Additionally, an [GitHub issue](https://help.github.com/en/enterprise/2.15/user/articles/creating-an-issue)
that outlines the feature that is to be implemented
or the bug to be fixed should be made. This issue will be closed in a 
[pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request) 
when the update is finished.

For other types of software changes (e.g. documentation) a similar structure may be used.

### Creating the branch
Before creating a branch, a 
[GitHub issue](https://help.github.com/en/enterprise/2.15/user/articles/creating-an-issue) 
should be made, as mentioned above.
The following command are used to create a new branch for a feature or bugfix:

```
git checkout develop
git checkout -b feature-demux_input_verify
```

From this point, the feature may be developed and tested as is seen fit.

# Merging the feature back into develop
Once the feature has been completed, the branch that was created 
earlier should be merged back into the develop branch. 
This is accomplished by submitting a [pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request) to merge the 
feature branch into PepSIRF's develop branch. [Zwfink](https://github.com/zwfink) 
should be requested as a reviewer.
Once approved, the feature may be merged.

**Note**: once a level of comfortability with PepSIRF development is reached, 
the branch may be merged manually with these commands:
```
git checkout develop
git merge --no-ff develop 
```

**It is important that you close the issue for this feature in this merge commit**.
See information about closing issues in commit messages [here](https://help.github.com/en/github/managing-your-work-on-github/linking-a-pull-request-to-an-issue).

After the branch that held the feature has been merged, it may be deleted.
After merging the branch, a note should be added to the [CHANGELOG.md](CHANGELOG.md) file under the
```Unreleased``` section giving a brief overview of the changes that were made. 

# Updating version number, creating a release
After one or more features have been implemented as described above, an official 
version may be released. The version number of the new release should be created as 
described in the semver specification.

To create a new release that is ready to be merged back into master, the following 
should be done: 

- Move the items that are in the ```Unreleased``` section of the CHANGELOG file into 
a new version, as is shown in the existing CHANGELOG file. The ```Unreleased``` heading
should remain, but will be empty.

- Update the current version listed in PepSIRF's main [README.md](README.md) file.

- Update the definition of ```PEPSIRF_VERSION``` located in 
[pepsirf_version.h](include/modules/core/pepsirf_version.h).


The release is now ready to be be merged into master.

# Merging back into master
Once a release is ready to be merged into master, the following steps may be done:

```
git checkout master
git merge master --no-ff 
git tag -a VERSION_NUMBER
git push 
git push --follow-tags
git checkout develop
git merge master
```

The new version of PepSIRF has been officially released! 
