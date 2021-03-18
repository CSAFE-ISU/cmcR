## Resubmission

This is a resubmission.

The previous version was rejected with the following comments:

"It seems we need to remind you of the CRAN policy:

'Packages which use Internet resources should fail gracefully with an informative message
if the resource is not available or has changed (and not give a check warning nor error).'

This needs correction whether or not the resource recovers."

Explanation of changes: Certain builds of the package (on r-devel-windows-ix86+x86_64 and r-release-windows-ix86+x86_64) produced a "cannot open URL" error due to a test that relies on downloading files from a U.S. government database. It's unclear whether this was due to a time-out, a lack of internet access, or some other reason. In any case, this (minor) test is now skipped on CRAN.

## Test environments
* local R installation: release
* ubuntu 16.04 (on travis-ci): release
* win-builder: devel
* GitHub Actions (windows): release
* GitHub Actions (ubuntu-20.04): release, devel

## R CMD check results

0 errors v | 1 warning | 0 notes v

Warning:

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1088/0957-0233/25/6/064005
    From: inst/doc/decisionRuleDescription.html
    Status: Error
    Message: libcurl error code 35:
      	schannel: next InitializeSecurityContext failed: SEC_E_ILLEGAL_MESSAGE (0x80090326) - This error usually occurs when a fatal SSL/TLS alert is received (e.g. handshake failed).
      	
Explanation of warning: This is not an invalid URL. I changed the URL used previously to link to this particular article to the more permanent DOI link, yet the warning persists.

#Previous cran-comments


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

