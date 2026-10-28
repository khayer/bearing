---
name: Test report
about: Report how installing and running BEARING on the practice data went (docs/TRY_BEARING.md)
title: "Test report: <your OS> - <PASS or FAIL>"
labels: test-report
---

<!--
  Issue template: test report for docs/TRY_BEARING.md
  Authors: Katharina E. Hayer with Claude (Anthropic), co-created
  Version: 1.0
  Date: 2026-10-28
-->

**OS / chip:** (for example, macOS 14 / Apple M2, Ubuntu 22.04, WSL2)

**Install route:** pip / conda

**Version line** (output of the python one-liner in TRY_BEARING.md step 5):

```
paste here
```

**pytest -q summary line:**

```
paste here
```

**Practice run (`bash workflow/smoke_test.sh --with-perm`):** PASS / FAIL

If FAIL, paste the last ~30 lines:

```
paste here
```

**Rough time per step** (install / tests / practice run):

**Anything confusing or missing in the instructions:**
