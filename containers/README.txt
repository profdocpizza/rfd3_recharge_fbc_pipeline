Container strategy
- One image per pipeline step/tool.
- Pin image by digest in production profiles.
- Keep a mutable :dev tag for local iteration.
- Rebuild only changed tool image when logic changes.
- Preserve old digests for historical run reproducibility.

Suggested image set
- hotspot_selector
- rfd3
- proteinrecharge
- freebindcraft
- mutant_zoo
- adapter/aggregate (python utility image)

