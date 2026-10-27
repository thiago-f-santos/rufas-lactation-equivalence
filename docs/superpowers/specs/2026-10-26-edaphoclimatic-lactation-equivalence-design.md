# Especificação de Design: Equivalência Edafoclimática da Curva de Lactação (EUA ↔ Brasil)

**Data:** 2026-10-26  
**Escopo:** Projeto autônomo `/home/thiago/Projects/rufas-lactation-equivalence/` e integração limpa em `/home/thiago/Projects/RuFaS/`

---

## 1. Objetivo e Contexto

No modelo RuFaS (*Ruminant Farm Systems*), a curva de lactação de Wood ($Y_t = l \cdot t^m \cdot e^{-nt}$) recebe ajustes regionais aditivos $(\Delta l, \Delta m, \Delta n)$ calibrados para 13 macrorregiões leiteiras dos Estados Unidos a partir dos códigos estaduais FIPS (`1` a `56`), conforme publicado por **Li, Reed & Cabrera (2022/2023)**.

Para simular rebanhos leiteiros nos 27 estados brasileiros (UFs do IBGE, códigos `11` a `53`) com fundamentação científica, transparência e reprodutibilidade:
1. Decompomos o perfil de clima e solo dos estados americanos mapeados no RuFaS e das 27 UFs brasileiras em **6 características biofísicas contínuas**.
2. Quantificamos e normalizamos cada característica em escala decimal adimensional $[0.0, 1.0]$.
3. Calculamos o **Índice de Similaridade Ambiental de Gower** ($S_{\text{Gower}} \in [0, 1]$) entre cada UF brasileira e os estados americanos, selecionando a região leiteira do estado americano mais próximo (**Top-1 Nearest Homologue**).
4. Executamos esse processamento de forma **atômica** em um projeto apartado (`rufas-lactation-equivalence/`), exportando apenas o mapeamento resultante (`brazil_state_to_region_mapping`) para o RuFaS e resolvendo a colisão numérica entre códigos IBGE de 2 dígitos e códigos FIPS estaduais dos EUA.

---

## 2. Fundamentação Científica e Modelo Matemático

### 2.1 Referências Bibliográficas Primárias
1. **Li, M., Reed, K. F., & Cabrera, V. E. (2022).** *Investigating the effect of temporal, geographic, and management factors on US Holstein lactation curve parameters.* Journal of Dairy Science, 105(9), 7525–7538. (Também *Li et al., 2023, J. Dairy Sci., 106(9), 6232–6248*).
2. **Bohmanova, J., Misztal, I., & Cole, J. B. (2007).** *Temperature-humidity indices as indicators of milk production losses due to heat stress.* Journal of Dairy Science, 90(4), 1947–1956. (Equação clássica NRC 1971 de THI).
3. **Fick, S. E., & Hijmans, R. J. (2017).** *WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas.* International Journal of Climatology, 37(12), 4302–4315.
4. **Fischer, G., et al. (2021).** *Global Agro-Ecological Zones (GAEZ v4) – Model Documentation.* FAO & IIASA.
5. **Poggio, L., et al. (2021).** *SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty.* SOIL, 7(1), 217–240.
6. **Saxton, K. E., & Rawls, W. J. (2006).** *Soil water characteristic estimates by texture and organic matter for hydrologic solutions.* Soil Science Society of America Journal, 70(5), 1569–1578.
7. **Gower, J. C. (1971).** *A general coefficient of similarity and some of its properties.* Biometrics, 27(4), 857–871.
8. **Ramirez-Villegas, J., Lau, C., Köhler, A. K., Signer, J., Jarvis, A., Arnell, N., Osborne, T., & Hooker, J. (2011).** *Climate analogues: finding tomorrow's agriculture today.* CGIAR Research Program on Climate Change, Agriculture and Food Security (CCAFS), Working Paper No. 12.

### 2.2 Equações Biofísicas Auxiliares

1. **Índice de Temperatura e Umidade (THI — NRC, 1971; Bohmanova et al., 2007):**
   Dado $T$ (temperatura média anual em °C) e $RH$ (umidade relativa média anual em %, $0 \le RH \le 100$):
   $$\text{THI} = (1.8 \cdot T + 32.0) - (0.55 - 0.0055 \cdot RH) \cdot (1.8 \cdot T - 26.0)$$

2. **Capacidade de Água Disponível no Solo (PAWC — Saxton & Rawls, 2006):**
   Dadas as frações texturais superficiais (0–30 cm) de areia $S \in [0, 1]$, argila $C \in [0, 1]$ e matéria orgânica $OM \in [0, 1]$ (onde $OM = \text{SOC}_{\text{g/kg}} \times 0.001724$):
   - Umidade no ponto de murcha permanente ($\theta_{1500}$, $-1500\text{ kPa}$):
     $$\theta_{1500t} = -0.024 S + 0.487 C + 0.006 OM + 0.005(S \cdot OM) - 0.013(C \cdot OM) + 0.068(S \cdot C) + 0.031$$
     $$\theta_{1500} = \theta_{1500t} + (0.14 \cdot \theta_{1500t} - 0.02)$$
   - Umidade na capacidade de campo ($\theta_{33}$, $-33\text{ kPa}$):
     $$\theta_{33t} = -0.251 S + 0.195 C + 0.011 OM + 0.006(S \cdot OM) - 0.027(C \cdot OM) + 0.452(S \cdot C) + 0.299$$
     $$\theta_{33} = \theta_{33t} + \left(1.283 \cdot \theta_{33t}^2 - 0.374 \cdot \theta_{33t} - 0.015\right)$$
   - Capacidade de Água Disponível ($\text{PAWC}$, $\text{cm}^3/\text{cm}^3$):
     $$\text{PAWC} = \max(0.01, \; \theta_{33} - \theta_{1500})$$

### 2.3 Decomposição nas 6 Características Normalizadas em $[0.0, 1.0]$

Cada região/estado $i$ é representado pelo vetor $\mathbf{x}_i = (x_{i,1}, x_{i,2}, x_{i,3}, x_{i,4}, x_{i,5}, x_{i,6}) \in [0.0, 1.0]^6$ usando transformação Min-Max truncada nos limites biofísicos globais $[B_{\min, k}, B_{\max, k}]$:

$$x_{i,k} = \text{clip}\left(\frac{X_{i,k} - B_{\min, k}}{B_{\max, k} - B_{\min, k}}, \; 0.0, \; 1.0\right)$$

| Índice | Chave no Código | Variável Física $X_{i,k}$ | Limite Mínimo $B_{\min, k}$ | Limite Máximo $B_{\max, k}$ |
| :--- | :--- | :--- | :--- | :--- |
| $x_1$ | `thi_index` | $\text{THI}$ médio anual (adimensional) | $45.0$ | $82.0$ |
| $x_2$ | `temp_seasonality` | $\Delta T = T_{\text{mês mais quente}} - T_{\text{mês mais frio}}$ (°C) | $1.0$ | $32.0$ |
| $x_3$ | `moisture_index` | Índice Hídrico $AI = P_{\text{anual}} / \text{PET}_{\text{anual}}$ | $0.10$ | $2.50$ |
| $x_4$ | `soil_clay_fraction` | Argila 0–30 cm (% em massa) | $5.0$ | $75.0$ |
| $x_5$ | `soil_pawc_index` | $\text{PAWC}$ de *Saxton & Rawls* ($\text{cm}^3/\text{cm}^3$) | $0.04$ | $0.22$ |
| $x_6$ | `soil_fertility_index` | $0.5 \cdot \text{norm}(\text{pH}, 4.2, 7.8) + 0.5 \cdot \text{norm}(\text{SOC}_{\text{g/kg}}, 5.0, 50.0)$ | $0.0$ | $1.0$ |

### 2.4 Métrica de Similaridade e Regra de Decisão Top-1

Para uma UF brasileira $b \in \{11, \dots, 53\}$ e um estado americano $u$ cujo mapeamento em `state_to_region_mapping` seja uma das 13 regiões leiteiras ativas ($R(u) \neq \text{"None"}$):

1. **Similaridade de Gower ($S_{\text{Gower}} \in [0, 1]$):**
   $$S_{\text{Gower}}(b, u) = 1.0 - \frac{1}{6} \sum_{k=1}^{6} |x_{b,k} - x_{u,k}|$$
2. **Similaridade Euclidiana Normalizada ($S_{\text{Eucl}} \in [0, 1]$, critério de auditoria e desempate):**
   $$S_{\text{Eucl}}(b, u) = 1.0 - \sqrt{\frac{1}{6} \sum_{k=1}^{6} (x_{b,k} - x_{u,k})^2}$$
3. **Seleção Top-1 Determinística:**
   Ordenam-se os estados americanos candidatos por $(S_{\text{Gower}}(b, u), S_{\text{Eucl}}(b, u), -\text{FIPS}_u)$ em ordem decrescente. O primeiro colocado $u^*$ define a região equivalente $R(u^*)$ atribuída à UF $b$.

---

## 3. Arquitetura do Projeto Apartado (`rufas-lactation-equivalence/`)

### 3.1 Coleta Automatizada de Dados Científicos (`fetcher.py`)
Para garantir transparência total e eliminar valores manuais arbitrários, os perfis brutos de clima e solo dos 48 estados americanos e das 27 UFs brasileiras são obtidos via código a partir de coordenadas geográficas de referência (`data/reference_coordinates.json`) consultando duas APIs científicas oficiais (sem necessidade de chave de API):
1. **Clima (NASA POWER Agroclimatology API — *Stackhouse et al., 2018*):**
   - Endpoint: `https://power.larc.nasa.gov/api/temporal/climatology/point?parameters=T2M,T2M_MAX,T2M_MIN,RH2M,PRECTOTCORR&community=AG`
   - Extrai: Temperatura média anual ($T_{\text{mean}}$, °C), Umidade relativa média anual ($RH_{\text{mean}}$, %), Sazonalidade térmica ($\Delta T = \max_m(T_{2M,m}) - \min_m(T_{2M,m})$, °C), Precipitação anual acumulada ($P_{\text{anual}} = \text{PRECTOTCORR}_{\text{ANN}} \times 365$, mm/ano) e Evapotranspiração Potencial anual ($\text{PET}_{\text{anual}}$, mm/ano).
2. **Solo 0–30 cm (ISRIC SoilGrids v2.0 OGC WCS 2.0.1 — *Poggio et al., 2021*):**
   - Endpoint: `https://maps.isric.org/mapserv?map=/map/<prop>.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage`
   - Propriedades consultadas na camada superficial (`0-5cm`, `5-15cm`, `15-30cm`, ponderadas por espessura $5:10:15$):
     - `clay` ($\text{g/kg} \div 10 \to \%$ em massa)
     - `sand` ($\text{g/kg} \div 10 \to \%$ em massa)
     - `phh2o` ($\text{pH} \times 10 \div 10 \to \text{pH em H}_2\text{O}$)
     - `soc` ($\text{dg/kg} \div 10 \to \text{g/kg de Carbono Orgânico}$)
   - Os coverages GeoTIFF 4×4 pixels em torno da coordenada são decodificados diretamente via biblioteca padrão do Python (`zlib` + `struct`), sem depender de GDAL/rasterio, e persistidos atomicamente em `data/us_states_edaphoclimatic.json` e `data/brazil_ufs_edaphoclimatic.json` com metadados de proveniência.

### 3.2 Estrutura de Arquivos
```text
/home/thiago/Projects/rufas-lactation-equivalence/
├── pyproject.toml
├── README.md
├── data/
│   ├── reference_coordinates.json
│   ├── us_states_edaphoclimatic.json
│   ├── brazil_ufs_edaphoclimatic.json
│   └── rufas_region_adjustments.json
├── src/
│   └── rufas_lactation_equivalence/
│       ├── __init__.py
│       ├── models.py
│       ├── biophysics.py
│       ├── fetcher.py
│       ├── similarity.py
│       ├── atomic_io.py
│       └── cli.py
├── output/
│   ├── equivalence_audit_report.json
│   └── equivalence_audit_report.md
└── tests/
    ├── test_biophysics.py
    ├── test_fetcher.py
    ├── test_similarity.py
    └── test_atomic_pipeline.py
```

### 3.3 Garantias de Execução Atômica (`atomic_io.py`)
Toda gravação em disco (geração dos relatórios em `output/` e exportação de `brazil_state_to_region_mapping` para `RuFaS/input/data/animal/lactation_curve_adjustment_inputs.json`) utiliza escrita transacional:
1. Serializa e valida todo o payload em memória antes de abrir qualquer descritor de arquivo.
2. Cria um `NamedTemporaryFile` no mesmo diretório do arquivo alvo (mesmo mount point / inode table).
3. Escreve os dados, executa `flush()` e `os.fsync()`.
4. Substitui o arquivo final via `os.replace(temp_path, target_path)`, garantindo atomicidade em nível de sistema operacional.

---

## 4. Contrato de Integração no `RuFaS/`

### 4.1 `input/data/animal/lactation_curve_adjustment_inputs.json`
- A chave `"state_to_region_mapping"` permanece inalterada (mapeando códigos FIPS dos EUA `"1"` a `"56"`).
- Adiciona-se a chave `"brazil_state_to_region_mapping"` mapeando os 27 códigos de UF do IBGE (`"11"`, `"12"`, `"13"`, `"14"`, `"15"`, `"16"`, `"17"`, `"21"`, `"22"`, `"23"`, `"24"`, `"25"`, `"26"`, `"27"`, `"28"`, `"29"`, `"31"`, `"32"`, `"33"`, `"35"`, `"41"`, `"42"`, `"43"`, `"50"`, `"51"`, `"52"`, `"53"`) para as respectivas regiões equivalentes do RuFaS calculadas pelo pipeline.

### 4.2 `RUFAS/biophysical/animal/milk/lactation_curve.py`
Em `LactationCurve.set_lactation_parameters`:
- Se `str(country).upper() == "BRA"`:
  Obtém `region_mapping = lactation_inputs.get("brazil_state_to_region_mapping", lactation_inputs.get("state_to_region_mapping", {}))`.
- Caso contrário (`"USA"`):
  Obtém `region_mapping = lactation_inputs.get("state_to_region_mapping", {})`.
- Isso garante que códigos IBGE como `31` (Minas Gerais) ou `3148004` (Patos de Minas/MG) consultem exclusivamente o mapa brasileiro, sem jamais colidir com o FIPS `31` (Nebraska) do mapa americano.

### 4.3 `RUFAS/input/metadata/properties/default.json`
- Adiciona-se a propriedade opcional `"brazil_state_to_region_mapping"` (`"type": "object"`, `"required": false`) dentro da seção `lactation_properties`, definindo os 27 códigos de UF do IBGE como `"type": "string"`.

---

## 5. Plano de Testes e Verificação

1. **Testes Unitários no `rufas-lactation-equivalence/`:**
   - Exatidão das fórmulas de THI e PAWC contra valores de referência da literatura.
   - Invariante de domínio: todas as 6 coordenadas de todos os estados americanos e brasileiros pertencem estritamente a $[0.0, 1.0]$.
   - Propriedades métricas de Gower ($S(a, a) = 1.0$, $S(a, b) = S(b, a)$, $0.0 \le S(a, b) \le 1.0$).
   - Atomicidade: simulação de falha de escrita não corrompe o arquivo alvo; execução com `--export-rufas-json` preserva todas as chaves existentes do JSON do RuFaS e atualiza apenas `brazil_state_to_region_mapping`.
2. **Testes Unitários e de Integração no `RuFaS/`:**
   - Validação de schema do `InputManager` carregando `lactation_curve_adjustment_inputs.json` atualizado.
   - Testes em `tests/test_biophysical/test_animal/test_milk_production/test_lactation_curve.py` verificando que `country="USA", region_code=31000` (Nebraska) e `country="BRA", region_code=3148004` (Minas Gerais) resolvem de forma independente a partir de seus respectivos dicionários de mapeamento.
