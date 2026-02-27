# Guidelines for creating meaningful issue reports

This section gives some started on how to generate useful issues on the
GROMACS [issue tracker](). The information
here comes to a large extent directly from there, to help you in
preparing your reports.

## What to report

Please only report issues you have confirmed to be caused by
GROMACS behaving in an unintended way, and
that you have investigated to the best of your ability. If you have
large simulations fail at some point, try to also trigger the problem
with smaller test cases that are more easily debuggable.

Bugs resulting from the use third-party software should be investigated
first to make sure that the fault is in
GROMACS and not in other parts of the
toolchain.

Please do not submit generic issues resulting from system instabilities
and systems `blowing-up`.

## What should be included

The report should include a general description of the problem with
GROMACS indicating both the expected behaviour
and the actual outcome. If the issue causes program crashes, the report
should indicate where the crash happens and if possible include the
stack trace right up to the crash.

All bugs should include the necessary information for the developers to
reproduce the errors, including if needed minimal input files (\*tpr,
\*top, \*mdp, etc), run commands or minimal version of run scripts, how
you compiled GROMACS and if possible the
system architecture.

The emphasis should be on having a *minimal* working example that is
easy to follow for the developers, that does not result in any warnings
or errors in itself. If your example generates errors, your issue will
likely not be considered as *real*, or at the minimum it will be much
harder to analyse to find the actual issue.

If your inputs are sensitive, then it is possible to create private
[issues](https://gitlab.com/gromacs/gromacs/-/issues/) so that the
developer team can have access to solve the problem, while preventing
widespread visibility on the internet.

## Supporting the developers

In general you should be able to answer questions posed to you by the
developers working on the program, if you want to help them in fixing
the bug you found. This may include things such as explaining run
scripts or simulation set-up, as well as confirming issues with
different versions of the program and different combinations of
supported libraries and compilers.

Please refrain from setting things such as target version or deciding on
unreasonable priorities. If you decide to fix the issue on your own,
please adhere to the other standards mentioned on the related pages
`code-formatting` and `code-commitstyle`.

> **See also:** `contribute`


## General issue workflow

The general issue workflow is shown in the figure below:

![Sample procedure pathway for reported issues.](redmine-states.png)

Project maintainers will apply [Status
labels](https://gitlab.com/gromacs/gromacs/-/labels?search=status) as
the issue is processed.

- [Status::Accepted](https://gitlab.com/gromacs/gromacs/-/issues?label_name%5B%5D=Status%3A%3AAccepted):
  Bug confirmed / Desirable feature.
- [Status::In
  Progress](https://gitlab.com/gromacs/gromacs/-/issues?label_name%5B%5D=Status%3A%3AIn+Progress):
  Assignee starts to work.
- [Status::Blocked](https://gitlab.com/gromacs/gromacs/-/issues?label_name%5B%5D=Status%3A%3ABlocked):
  Progress requires feedback or other action.
- [Status::Rejected](https://gitlab.com/gromacs/gromacs/-/issues?label_name%5B%5D=Status%3A%3ARejected):
  Invalid report or not a desirable feature.
- [Status::Fix
  uploaded](https://gitlab.com/gromacs/gromacs/-/issues?label_name%5B%5D=Status%3A%3AFix+uploaded):
  Merge request is available for review
- [Status::Feedback-wanted](https://gitlab.com/gromacs/gromacs/-/issues?label_name%5B%5D=Status%3A%3AFeedback-wanted):
  Resolution pending additional feedback or response
- [Status::Resolved](https://gitlab.com/gromacs/gromacs/-/issues?label_name%5B%5D=Status%3A%3AResolved):
  The issue will be closed if there is no further discussion.
