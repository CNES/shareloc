# **Shareloc** **Contributing guide**.

1. [Bug report](#bug-report)
2. [Contributing workflow](#contributing-workflow)
3. [Contribution license agreement](#contribution-license-agreement)
4. [Coding guide](#coding-guide)
5. [Pre-commit validation](#pre-commit-validation)

# Bug report

Any proven or suspected malfunction should be traced in a bug report, the latter being an issue in the Shareloc repository.

**Don't hesitate to do so: It is best to open a bug report and quickly resolve it than to let a problem remains in the project.**
**Notifying the potential bugs is the first way for contributing to a software.**

In the problem description, be as accurate as possible. Include:
* The procedure used to initialize the environment
* The incriminated python function or class
* The content of the input and output data

# Contributing workflow

Any code modification requires a Pull/Merge Request. It is forbidden to push patches directly into master (this branch is protected).

It is recommended to open your Merge Request as soon as possible in order to inform the developers of your ongoing work.
Please add `WIP:` before your Merge Request title if your work is in progress: This prevents an accidental merge and informs the other developers of the unfinished state of your work.

The Merge Request shall have a short description of the proposed changes. If it is relative to an issue, you can signal it by adding `Closes xx` where xx is the reference number of the issue.

Likewise, if you work on a branch (which is recommended), prefix the branch's name by `xx-` in order to link it to the xx issue.

Shareloc Classical workflow is :
* Check Licence and sign [Contributor Licence Agrement](#contribution-license-agreement) (Individual or Corporate)
* Create an issue (or begin from an existing one)
* Create a Merge Request from the issue: a MR is created accordingly with "WIP:", "Closes xx" and associated "xx-name-issue" branch
* Modify Shareloc code from a local working directory or from the forge (less possibilities)
* Git add, commit and push from local working clone directory or from the forge directly
* Follow [Conventional commits](https://www.conventionalcommits.org/) specifications for commit messages
* Beware that pre-commit hooks can be installed for code analysis (see below pre-commit validation).
* Launch the tests with [pytest](https://pytest.org) on your modifications (or don't forget to add ones).
* When finished, change your Merge Request name (erase "WIP:" in title ) and ask to review the code.
 

# Contribution license agreement

Shareloc requires that contributors sign out a [Contributor License
Agreement](https://en.wikipedia.org/wiki/Contributor_License_Agreement). The
purpose of this CLA is to ensure that the project has the necessary ownership or
grants of rights over all contributions to allow them to distribute under the
chosen license (Apache License Version 2.0)

To accept your contribution, we need you to complete, sign and email to *cars [at]
cnes [dot] fr* an [Individual Contributor Licensing
Agreement](./docs/source/CLA/ICLA-Shareloc.doc) (ICLA) form and a
[Corporate Contributor Licensing
Agreement](./docs/source/CLA/CCLA-Shareloc.doc) (CCLA) form if you are
contributing on behalf of your company or another entity which retains copyright
for your contribution.

The copyright owner (or owner's agent) must be mentioned in headers of all
modified source files and also added to the [NOTICE
file](./NOTICE).


