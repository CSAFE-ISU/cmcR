## Resubmission

Dr. Ligges requested re-writing a DOI address in the Description file.

### Explanation of changes

I have rewritten the DOI address from "<doi:10.6028> (10.6028/jres.120.008)" to "<doi:10.6028/jres.120.008>"

## Version 0.1.10 Test environments
* local R installation: release
* win-builder: devel
* GitHub Actions (windows): release
* GitHub Actions (ubuntu-20.04): release, devel

#Previous cran-comments

## Resubmission

The cmcR tests failed due to recent changes to the tidyselect package.

### Explanation of changes

An author of tidyselect, Hadley Wickham, made the necessary changes to work with tidyselect.
See the following github commit for proof: https://github.com/CSAFE-ISU/cmcR/commit/974f5b9fb3ceea3956019ca9bc43d184f0ff4f1c
Tests are passing without notes, warnings, or errors on all test environments

## Version 0.1.10 Test environments
* local R installation: release
* win-builder: devel
* GitHub Actions (windows): release
* GitHub Actions (ubuntu-20.04): release, devel

## Resubmission

There was an "Additional issue" while running an example that included a donttest tag.
This issue has been resolved by changing this example tag to dontrun.

### Explanation of changes

Two links on the README.md file required updating. Specifically, the codecov address is now https://app.codecov.io/gh/CSAFE-ISU/cmcR?branch=master while the Travis-CI link has been removed.

The other URLs included in the note link to U.S. government websites. Each of these have been checked and link to the appropriate webpage.

## Version 0.1.8 Test environments
* local R installation: release
* ubuntu 16.04 (on travis-ci): release
* win-builder: devel
* GitHub Actions (windows): release
* GitHub Actions (ubuntu-20.04): release, devel

## Resubmission

The previous submission had a note that of (possibly) invalid URLs. As explained below, two of these thinks indeed needed updating while the rest have been noted previously and deemed valid.

### Explanation of changes

Two links on the README.md file required updating. Specifically, the codecov address is now https://app.codecov.io/gh/CSAFE-ISU/cmcR?branch=master while the Travis-CI link has been removed.

The other URLs included in the note link to U.S. government websites. Each of these have been checked and link to the appropriate webpage.

## Version 0.1.8 Test environments
* local R installation: release
* ubuntu 16.04 (on travis-ci): release
* win-builder: devel
* GitHub Actions (windows): release
* GitHub Actions (ubuntu-20.04): release, devel

## Resubmission

The previous submission did not automatically pass the CRAN checks because I had not updated the package version. The package was previously archived on 2021-07-12 because it requires a previously archived package x3ptools that has since been reinstated to CRAN. No other changes have been made to the package other than updating the package version to 0.1.7.

## Version 0.1.7 Test environments
* local R installation: release
* ubuntu 16.04 (on travis-ci): release
* win-builder: devel
* GitHub Actions (windows): release
* GitHub Actions (ubuntu-20.04): release, devel

## Version 0.1.7 R CMD check results

0 errors v | 0 warnings | 0 notes

## Resubmission

This is a resubmission of cmcR version 0.1.6. cmcR was removed from CRAN because a package upon which it relies, x3ptools, was removed from CRAN on 2021-06-27. The x3ptools issue has been fixed and x3ptools is now available on CRAN. There have been no changes to cmcR since its initial acceptance onto CRAN on 2021-04-05, so this submission is merely meant to reinstate it on CRAN.

## Version 0.1.6 Test environments
* local R installation: release
* ubuntu 16.04 (on travis-ci): release
* win-builder: devel
* GitHub Actions (windows): release
* GitHub Actions (ubuntu-20.04): release, devel

## Version 0.1.6 R CMD check results

0 errors v | 0 warnings | 0 notes

## Resubmission

This is a resubmission following a failed Version 0.1.5 resubmission. Version 0.1.5 was an attempt to fix a problem identified in an email from Professor Ripley regarding a failing test in version 0.1.4 of cmcR (see details below). Version 0.1.5 fixed the problem identified by Professor Ripley, yet did not pass the automatic CRAN tests 9see details below). Version 0.1.6 addresses the causes (2 Notes) for the failed Version 0.1.5 submission.

## Version 0.1.6 Test environments
* local R installation: release
* ubuntu 16.04 (on travis-ci): release
* win-builder: devel
* GitHub Actions (windows): release
* GitHub Actions (ubuntu-20.04): release, devel

## Version 0.1.6 R CMD check results

0 errors v | 0 warnings | 0 notes

### Version 0.1.5 submission 2 Notes

Following is an explanation of the 2 Notes from the attempted submission of version 0.1.5.

#### Note 1

"Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1088/0957-0233/25/6/064005
    From: inst/doc/decisionRuleDescription.html
    Status: 403
    Message: Forbidden
  URL: https://nvlpubs.nist.gov/nistpubs/jres/120/jres.120.008.pdf
    From: inst/doc/decisionRuleDescription.html
    Status: 404
    Message: Not Found
  URL: https://tsapps.nist.gov/NRBTD/
    From: inst/doc/cmcR_plotReproduction.html
    Status: 404
    Message: Not Found
  URL: https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/Details/2d9cc51f-6f66-40a0-973a-a9292dbee36d
    From: man/fadulData_processed.Rd
    Status: 404
    Message: Not Found
  URL: https://tsapps.nist.gov/NRBTD/Studies/Firearm/Details/12e0f761-2528-4e7b-8600-360bbb788537
    From: inst/doc/cmcR_plotReproduction.html
    Status: 404
    Message: Not Found
  URL: https://tsapps.nist.gov/NRBTD/Studies/Firearm/Details/681f3cdf-8b1c-418a-af71-f52fd235e3da
    From: inst/doc/cmcR_plotReproduction.html
    Status: 404
    Message: Not Found
  URL: https://tsapps.nist.gov/NRBTD/Studies/Search
    From: README.md
    Status: 404
    Message: Not Found
  URL: https://tsapps.nist.gov/NRBTD/Studies/Studies/Details/e7a8aab8-8d5a-44ac-b2be-f0de7c2ca505?nm=True&mt=1&m=3&sp=1
    From: inst/doc/decisionRuleDescription.html
    Status: 404
    Message: Not Found
  URL: https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=910868
    From: inst/doc/decisionRuleDescription.html
    Status: 404
    Message: Not Found
  URL: https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=911193
    From: DESCRIPTION
          man/decision_CMC.Rd
          README.md
    Status: 404
    Message: Not Found"

Explanation of Note 1: None of these are invalid URLs. The one DOI is a valid DOI while the other URLs link to a U.S. government agency database.

#### Note 2

"* checking top-level files ... NOTE
Non-standard file/directory found at top level:
  ‘tests_NOT_CRAN’"
  
Explanation of Note 2: This was a non-standard folder structure used to avoid running a test that failed on CRAN (see comment from Professor Ripley concerning version 0.1.4 submission below). I realize that this was not an effective solution to skip this test, and have since implemented a solution using testthat::skip_on_cran.

### Version 0.1.4 comments from Professor Ripley

A previous version, version 0.1.4, was rejected with the following comments from Professor Ripley:

"It seems we need to remind you of the CRAN policy:

'Packages which use Internet resources should fail gracefully with an informative message
if the resource is not available or has changed (and not give a check warning nor error).'

This needs correction whether or not the resource recovers."

Explanation of changes addressing Professor Ripley's comments: Certain builds of the package (on r-devel-windows-ix86+x86_64 and r-release-windows-ix86+x86_64) produced a "cannot open URL" error due to a test that relies on downloading files from a U.S. government database. It's unclear from the log files specifically why this error occurs and we have been unable to replicate the problem in any other test environment. This test is now skipped on CRAN.

## Resubmission

This is a minor, yet necessary, package update that fixes bugs in the package's plotting functions x3pListPlot and cmcPlot.

## Test environments
* local R installation: release
* ubuntu 16.04 (on travis-ci): release
* win-builder: devel
* GitHub Actions (windows): release
* GitHub Actions (ubuntu-20.04): release, devel

## R CMD check results

0 errors v | 0 warnings v | 0 notes v

R CMC check succeeded

## Resubmission
This is a resubmission.
The previous submission was rejected with the following comments:

"Please reduce the length of the title to less than 65 characters."

* The title is now 58 characters long.

"Please always write package names, software names and API (application
programming interface) names in single quotes in title and description.
e.g: --> 'cmcR'"

* I have added quotations around 'Congruent Matching Cells.' To my knowledge, there is not another package name, software name, or API referenced in the title or description.

"Please write references in the description of the DESCRIPTION file in
the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: authors (year) <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")"

* I have added to the description appropriate links to Song (2013) (note that the URL is the only available identifying information, so I had to use it instead of a DOI, etc.) and Tong et al. (2015) (which does have a DOI).

"\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user.
Does not seem necessary.

Please unwrap the examples if they are executable in < 5 sec, or create
additionally small toy examples to allow automatic testing.
(You could also replace \dontrun{} with \donttest, if it takes longer
than 5 sec to be executed, but it would be preferable to have automatic
checks for functions. Otherwise, you can also write some tests.)"
 
* I have removed \dontrun{} from examples that take < 5 seconds to run (total example run time on my local machine is now 46.5 seconds). The remaining examples involve downloading large 3D scans from a remote database and processing them and thus take more than 5 seconds. As this package is intended to be used by forensic practitioners (who may not have a thorough understanding of navigating R documentation files), my intention is to make each of these time-intensive examples "self-contained." These scans are too large in their unprocessed format to include with the package (at least to keep the package size below 5 MB), but I have included 2 processed scans in the data directory that are of reasonable size. The tests I have in the tests directory perform the automatic testing that Dr. Seyer references in his comment, so removing \dontrun{} around the time-intensive examples for the sake of testing would be redundant.


## Test environments
* local R installation: release
* ubuntu 16.04 (on travis-ci): release
* win-builder: devel
* GitHub Actions (windows): release
* GitHub Actions (ubuntu-20.04): release, devel

## R CMD check results

0 errors | 0 warnings | 0 notes

#Previous cran-comments

## Resubmission
This is a resubmission.
The previous submission was rejected with the following message from Professor Ligges:

"We see

* checking re-building of vignette outputs ... [454s] OK

which is somewhat long as the test also run 400 sec on a bi-arch platform.


Please reduce."

The build time for the vignettes has been reduced by changing the way certain plots are constructed (which was the major bottleneck in the previous version of the package).

## Test environments
* local R installation: release
* ubuntu 16.04 (on travis-ci): release
* win-builder: devel
* GitHub Actions (windows): release
* GitHub Actions (ubuntu-20.04): release, devel

## R CMD check results

0 errors | 0 warnings | 1 note

1 note is:
Possibly mis-spelled words in DESCRIPTION:
  al (13:55)
  cmcR (16:9)
  et (13:52)
  pre (14:30)

This is the same note as in the previous submission. The short of it is that these are not mis-spelled words.

## Resubmission
This is a resubmission.
The previous submission failed the automatic check due to a warning thrown while running the vignette.
The cause of the warnings has been fixed.

## Test environments
* local R installation: release
* ubuntu 16.04 (on travis-ci): release
* win-builder: devel
* GitHub Actions (windows): release
* GitHub Actions (ubuntu-20.04): release, devel

## R CMD check results

0 errors | 0 warnings | 1 note

1 note is:
Possibly mis-spelled words in DESCRIPTION:
  al (13:55)
  cmcR (16:9)
  et (13:52)
  pre (14:30)

Justification for the note: the usage of these is intentional.
* "et" and "al" are used to refer to "and others" since the package contains an implementation of a method from Tong "et al." (2015), which needs to be mentioned in the description as the original paper is not my own
* "cmcR" is the name of the package
* "pre" is because the package contains various "pre, inter, and post-processing" functions

# Previous cran-comments
## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel)
* ubuntu 20.04 (on GitHub Actions), R 4.0.3
* windows 10.0.17763 (on GitHub Actions), R 4.0.3

## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel)
* ubuntu 20.04 (on GitHub Actions), R 4.0.3
* windows 10.0.17763 (on GitHub Actions), R 4.0.3

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

