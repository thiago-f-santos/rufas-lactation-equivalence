# RuFaS Lactation Curve Edaphoclimatic Equivalence (USA FIPS ↔ Brazil IBGE)

Pacote científico determinístico e de execução atômica para decompor perfis de clima e solo dos estados americanos (códigos FIPS) e das 27 Unidades Federativas do Brasil (códigos IBGE) em vetores normalizados $[0.0, 1.0]^6$, identificando a região leiteira homóloga mais próxima nos EUA (**Top-1 Nearest Homologue**) via Coeficiente de Similaridade de Gower (1971) para parametrizar a curva de lactação de Wood no **RuFaS**.

## Execução Atômica

```bash
# Gerar relatórios científicos de auditoria (JSON + Markdown) e exportar para o RuFaS:
PYTHONPATH=src python3 -m rufas_lactation_equivalence.cli \
  --export-rufas-json ../RuFaS/input/data/animal/lactation_curve_adjustment_inputs.json
```

## Referências Científicas

1. **Li, M., Reed, K. F., & Cabrera, V. E. (2022).** *Investigating the effect of temporal, geographic, and management factors on US Holstein lactation curve parameters.* Journal of Dairy Science, 105(9), 7525–7538.
2. **Bohmanova, J., Misztal, I., & Cole, J. B. (2007).** *Temperature-humidity indices as indicators of milk production losses due to heat stress.* Journal of Dairy Science, 90(4), 1947–1956.
3. **Saxton, K. E., & Rawls, W. J. (2006).** *Soil water characteristic estimates by texture and organic matter for hydrologic solutions.* Soil Science Society of America Journal, 70(5), 1569–1578.
4. **Poggio, L., et al. (2021).** *SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty.* SOIL, 7(1), 217–240.
5. **Fick, S. E., & Hijmans, R. J. (2017).** *WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas.* International Journal of Climatology, 37(12), 4302–4315.
6. **Gower, J. C. (1971).** *A general coefficient of similarity and some of its properties.* Biometrics, 27(4), 857–871.
7. **Ramirez-Villegas, J., et al. (2011).** *Climate analogues: finding tomorrow's agriculture today.* CGIAR CCAFS Working Paper No. 12.
