Outside of the loop initialize 

```
L2L3Residual * L2L3 = new L2L3Residual(radius, etacut, dopPb);
```

Possibilities are radius = 3, etacut = 3 or 4 at the moment.

Then for each jet
```
get_corrected_pt(jtpt, jteta);
```