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

#Previous cran-comments
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

