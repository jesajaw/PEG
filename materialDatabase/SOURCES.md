# Material database

Refractive index tables used by PEG. One file per material, named `<Chemical Formular>.idx` (for example `Au.idx`, `SiO2.idx`).

## Where the data comes from

Henke tables from CXRO (Center for X-ray Optics, Lawrence Berkeley National Laboratory):

* Website: <https://henke.lbl.gov/optical_constants/>
* Index of refraction ($$\delta$$ and $$\beta$$): <https://henke.lbl.gov/optical_constants/getdb.html>
* Atomic scattering factors $$f1$$ and $$f2$$ per element: <https://henke.lbl.gov/optical_constants/asf.html>

Reference: B. L. Henke, E. M. Gullikson and J. C. Davis, *X-ray interactions: photoabsorption, scattering, transmission, and reflection at E = 50-30000 eV, Z = 1-92*, Atomic Data and Nuclear Data Tables **54**(2), 181-342 (1993).

## File format

```text
wavelength  delta   beta
```
Basically just three numbers per line, separated by blanks or tabs. The Refractive index notation:
$$ n = 1 - \delta + \mathrm{i}\beta, \beta \ge 0$$
No headers or comments, waveölength is in `nm` and increasing/ascending and in range of 0.0413nm to 40nm (30eV to 30keV).
If you want to use electronvolts instead, PEG uses:
$$ E = \frac{\mathrm{h}c}{\lambda} \approx\frac{1239.84172}{\lambda} $$
to convert it.