# Especificação de Design: Dashboard Estático de Equivalência Edafoclimática e Curvas de Lactação (Brasil ↔ EUA)

**Data:** 2026-10-27  
**Escopo:** Módulo visual e relatório interativo em `/home/thiago/Projects/rufas-lactation-equivalence/docs/index.html`  
**Status:** Aprovado para Implementação  

---

## 1. Visão Geral e Contexto Científico

No projeto `rufas-lactation-equivalence`, consolidamos a conversão determinística das curvas de lactação do modelo **RuFaS** (*Ruminant Farm Systems*), parametrizadas originalmente para 13 macrorregiões leiteiras dos Estados Unidos a partir de códigos estaduais FIPS (*Li, Reed & Cabrera, 2022/2023*), para as 27 Unidades Federativas do Brasil (códigos IBGE `11` a `53`).

Essa conversão decompôs cada estado em um vetor normalizado $[0.0, 1.0]^6$ baseado em 6 características biofísicas fundamentais de solo e clima:
1. `thi_index`: Índice de Temperatura e Umidade (NRC 1971 / Bohmanova et al., 2007).
2. `temp_seasonality`: Sazonalidade Térmica ($\Delta T$, amplitude térmica anual).
3. `moisture_index`: Índice Hídrico ($P / \text{PET}$, balanço precipitação / evapotranspiração).
4. `soil_clay_fraction`: Teor de argila superficial 0–30 cm (%).
5. `soil_pawc_index`: Capacidade de retenção de água no solo (*Saxton & Rawls, 2006*).
6. `soil_fertility_index`: Fertilidade química integrada (pH em água + Carbono Orgânico do Solo).

Através do Coeficiente de Similaridade de Gower ($S_{\text{Gower}} \in [0, 1]$), cada UF brasileira foi associada ao seu **Top-1 Nearest Homologue** americano, herdando a região leiteira correspondente do RuFaS e seus ajustes aditivos $(\Delta l, \Delta m, \Delta n)$ na curva de lactação de Wood:
$$Y(t) = l \cdot t^m \cdot e^{-n \cdot t} \quad (t \in [1, 305] \text{ DIM})$$

O objetivo desta especificação é definir o design de uma página estática interativa e autocontida (`docs/index.html`), equipada com um mapa vetorial SVG do Brasil, capaz de exibir de forma clara e elegante, no clique e no hover:
- O estado brasileiro selecionado e sua localização geográfica.
- As características biofísicas de clima e solo que decompusemos (valores brutos reais e normalizados).
- Para qual clima/estado americano foi roteado, sua região no RuFaS e grau de similaridade.
- A curva de lactação específica modelada para a UF comparada com a curva base neutra do RuFaS.

---

## 2. Requisitos Funcionais e Experiência do Usuário

### 2.1 Mapa Interativo do Brasil em SVG
- **Representação Geográfica:** 27 Unidades Federativas representadas por polígonos vetoriais `<path>` derivados da malha geográfica oficial do IBGE.
- **Coloração por Região RuFaS:** Cada estado no mapa é colorido de acordo com sua região RuFaS equivalente:
  - **Delta:** Verde esmeralda (`#10b981`) — ex: MG, MT, GO, RJ, TO, AC, AM, AP, DF.
  - **Southeast:** Azul oceano (`#3b82f6`) — ex: SP, ES, RO, RR, PA, MA, PI, CE, PB, PE, MS.
  - **Mountain:** Âmbar / Laranja dourado (`#f59e0b`) — ex: RN, SE.
  - **Southern Plains:** Vermelho carmim (`#ef4444`) — ex: BA, AL.
  - **Appalachian:** Roxo violeta (`#8b5cf6`) — ex: PR, SC, RS.
- **Interação no Hover:**
  - O estado ganha brilho/destaque de borda.
  - Tooltip flutuante suave exibe:
    - UF e Nome do Estado.
    - Estado Homólogo dos EUA e sigla (ex: *Louisiana, LA*).
    - Região RuFaS atribuída e Similaridade de Gower ($S_{\text{Gower}}$).
    - Pico diário estimado (kg/dia e dia do pico).
- **Interação no Clique:**
  - Fixa o estado como selecionado (borda com contorno em destaque e sombra).
  - Atualiza todos os painéis de dados, gráficos e cards de métricas.
  - Sincroniza o seletor/dropdown de estados.

### 2.2 Curva de Lactação de Wood ($Y(t)$)
- **Cálculo da Curva:** Produção diária de leite para cada dia em lactação $t \in [1, 305]$:
  - $l = 19.90 + \Delta l$
  - $m = 0.247 + (\Delta m \times 10^{-2})$
  - $n = 0.003376 + (\Delta n \times 10^{-4})$
  - $Y(t) = l \cdot t^m \cdot e^{-n \cdot t}$
- **Gráfico Comparativo:**
  - Linha sólida colorida com a cor da região RuFaS para a UF selecionada.
  - Linha pontilhada neutra (cinza) para a curva padrão média dos EUA/RuFaS ($l=19.9, m=0.247, n=0.003376$).
  - Ponto de destaque no vértice da curva indicando o dia do pico ($t_{\text{peak}} = m/n$) e produção de pico ($Y_{\text{peak}} = l \cdot (m/n)^m \cdot e^{-m}$).
- **Indicadores Numéricos:**
  - Produção de Pico ($Y_{\text{peak}}$, kg/dia).
  - Dia do Pico ($t_{\text{peak}}$, DIM).
  - Produção Acumulada em 305 dias ($\sum_{t=1}^{305} Y(t)$, kg totais).
  - Parâmetros ajustados exatos ($l, m, n$) e seus deltas $(\Delta l, \Delta m, \Delta n)$.

### 2.3 Decomposição Edafoclimática dos 6 Fatores
- **Gráfico Radar / Barras Biofísicas:**
  - Mostra o vetor normalizado $[0.0, 1.0]$ da UF comparado ao vetor do estado americano homólogo.
- **Detalhamento das Grandezas Físicas Reais:**
  - Tabela/grid de cartões com os valores físicos reais:
    - Temperatura média anual (°C) e Sazonalidade ($\Delta T$, °C).
    - Umidade Relativa média anual (%) e Índice THI.
    - Precipitação anual (mm/ano) e Índice Hídrico ($P / \text{PET}$).
    - Argila e Areia na camada 0–30 cm (% em massa).
    - PAWC ($\text{cm}^3/\text{cm}^3$) de *Saxton & Rawls*.
    - pH em água e Carbono Orgânico do Solo (g/kg).

### 2.4 Tabela de Auditoria das 27 UFs
- Tabela localizada na parte inferior com todas as 27 UFs.
- Colunas: Código IBGE, UF, Estado, Região BR, Estado Homólogo EUA, Região RuFaS, Similaridade Gower, Similaridade Euclidiana, $\Delta l$, $\Delta m$, $\Delta n$, Pico (kg/d), Total 305d (kg).
- Campo de filtro por texto em tempo real e ordenação por clique no cabeçalho das colunas.
- Clicar em qualquer linha seleciona a UF no mapa e carrega os gráficos.

---

## 3. Arquitetura de Implementação e Código

### 3.1 Estrutura de Arquivos
```text
rufas-lactation-equivalence/
├── docs/
│   ├── index.html                               # Dashboard estático autocontido
│   └── superpowers/specs/2026-10-27-...md       # Este documento de especificação
├── scripts/
│   └── build_dashboard.py                       # Script gerador/atualizador do dashboard
├── tests/
│   └── test_dashboard_integrity.py              # Teste unitário de integridade do dashboard
├── data/
│   ├── brazil_ufs_edaphoclimatic.json
│   ├── us_states_edaphoclimatic.json
│   └── rufas_region_adjustments.json
└── output/
    ├── equivalence_audit_report.json
    └── equivalence_audit_report.md
```

### 3.2 Autocontenção e Portabilidade de `docs/index.html`
1. **Zero Dependência Local:** Todo o CSS é embutido na tag `<style>`. O layout é responsivo, construído com CSS Grid e Flexbox, funcionando perfeitamente em telas desktop e tablets/smartphones.
2. **Dados Inlined:** O payload estruturado em JSON com as 27 UFs, suas geometrias SVG, vetores edafoclimáticos e curvas de Wood é injetado diretamente no documento através de `<script id="dashboard-data" type="application/json">`.
3. **Compatibilidade de Execução:** Pode ser aberto diretamente via `file:///.../docs/index.html` com duplo clique no navegador, ou servido via GitHub Pages ou qualquer servidor web estático sem problemas de CORS.
4. **Gráficos Fluidos:** Carregamento do Chart.js via CDN (`https://cdn.jsdelivr.net/npm/chart.js`) para visualização responsiva de alta performance com tooltips nativos, com suporte a fallback visual nativo.

### 3.3 Script Construtor Determinístico (`scripts/build_dashboard.py`)
Para manter a sincronização perfeita entre os relatórios do pipeline e o arquivo HTML:
- O script lê `output/equivalence_audit_report.json` e `data/brazil_ufs_edaphoclimatic.json`.
- Obtém as geometrias oficiais simplificadas das 27 UFs do IBGE e projeta em coordenadas de tela SVG ($600 \times 600$ com projeção Mercator ajustada aos limites do Brasil: longitude $[-74.0, -34.8]$, latitude $[-33.8, +5.3]$).
- Computa os pontos da curva de Wood para cada uma das 27 UFs.
- Monta e salva `docs/index.html` de forma atômica.

---

## 4. Plano de Verificação e Testes

O teste automatizado `tests/test_dashboard_integrity.py` garantirá que:
1. O arquivo `docs/index.html` existe e seu tamanho é coerente (> 30 KB).
2. O mapa SVG contém exatamente 27 elementos `<path>` correspondentes às 27 UFs brasileiras.
3. Todas as 27 UFs possuem:
   - Código IBGE e sigla válidos.
   - Estado americano homólogo preenchido.
   - Ajustes de Wood $(\Delta l, \Delta m, \Delta n)$ condizentes com a região RuFaS.
   - Similaridade de Gower no intervalo $[0.75, 1.0]$.
4. Os cálculos matemáticos do dia de pico e produção total em 305 dias no HTML batem com as equações de *Li et al. (2022)*.
