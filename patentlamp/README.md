# PatentLAMP 🛡️🧬
**Autonomous Patent-Ready Molecular Engineering & Intellectual Property Engine for Isothermal Amplification Primers**

---

## 📌 Visão Geral & Ideia Central

O **PatentLAMP** é uma plataforma autônoma de bioinformática e engenharia molecular projetada para resolver a lacuna crítica entre o **desenho bioinformático de primers** (ex: via *NextLAMP*, *Primer3* ou *GLAPD*) e a sua **efetiva proteção patentária e transferência de tecnologia**.

Na maioria das legislações internacionais — e especialmente no Brasil sob a **Lei de Propriedade Industrial (LPI nº 9.279/1996, Art. 10, IX)** — sequências de DNA/RNA puras ou isoladas exatamente como existem na natureza **não são patenteáveis**. 

O **PatentLAMP** automatiza a transformação de primers "zerados" (naturais) em **moléculas quiméricas sintéticas não-naturais**, calculando e provando matematicamente o seu **Efeito Técnico Inesperado (LPI Art. 13)** em relação ao controle natural não-modificado e compilando todo o pacote de depósito legal no padrão **WIPO ST.26 XML**.

---

## 🎯 Objetivos Principais & Funcionalidades Implementadas

1. **Superar a Vedação de Patenteabilidade de DNA Natural (LPI Art. 10):**
   * Engenharia algorítmica de junções sintéticas não-naturais (*Self-Folding Linkers*, caudas sintéticas, ligantes não-nucleosídicos como HEG/Spacer C3 e modificações químicas LNA) nas fitas FIP e BIP.
   * Suporte universal a ligantes **nucleosídicos (DNA sintético A,T,C,G)** e **não-nucleosídicos químicos (`HEG`, `SPACER_C3`)**, com representação oficial `n` e anotações de recurso `misc_feature` no XML ST.26.
   * Suporte a ligantes customizados de qualquer sequência digitados diretamente pela linha de comando.

2. **Comprovar a Atividade Inventiva / Efeito Técnico Surpreendente (LPI Art. 13 & USPTO Non-Obviousness):**
   * **Ensaio Comparativo Pares In Silico:** Comparação pareada e direta entre o *Primer Sintético Modificado* e o *Primer Natural Controle* (antes da modificação).
   * **Termodinâmica Rígida de SantaLucia (1998):** Cálculo de energia livre ($\Delta G$) por matriz de alinhamento com deslocamento total (*full-length offset matrix alignment scanning*) a 65°C.
   * **Fator de Supressão Termodinâmica de Boltzmann:** Cálculo físico da constante $S_{\text{factor}} = e^{\Delta \Delta G / RT}$ sem multiplicadores estáticos ou números fictícios.
   * **Justificativas Condicionais Sem Contradição:** O texto de prova adapta-se matematicamente: se $\Delta \Delta G > 0$, destaca a supressão de auto-dimerização (Art. 13); se $\Delta \Delta G \le 0$, fundamenta-se estritamente na Não-Naturalidade Quimérica Sintética (Art. 10, IX).

3. **Automação do Pacote de Depósito de Patente Internacional (WIPO ST.26 & INPI):**
   * **WIPO Standard ST.26 v1.3 XML:** Geração de arquivo XML com validação de DTD estrita (`<INSDSeq_feature-table>`), anotações de `source` (`other DNA`, `synthetic construct`), `misc_feature` para ligantes e entrada própria de `SEQ ID NO` para o ligante sintético isolado.
   * **Minuta de Patente Markdown (`INPI_Patent_Specification_Draft.md`):** Redação em português técnico no padrão da Resolução INPI nº 248/2019, citando nominalmente cada `SEQ ID NO` no Quadro Reivindicatório (Produtos Quiméricos, Kit Diagnóstico e Método *In Vitro*).
   * Citação dinâmica de primers de alça (`LoopF` e `LoopB`), adaptando o texto de acordo com a composição do conjunto.

4. **Robustez de Software & Tratamento de Dados:**
   * Suporte ao processamento simultâneo de **Múltiplos Conjuntos de Primers (`Multi-Set Support`)** na CLI e no relatório descritivo.
   * Resolução de ligantes e chaves de dicionário **insensíveis a maiúsculas/minúsculas (*case-insensitive*)** com mapeamento de sinônimos (`spacer18`, `c3`, `heg`, `taaa`).
   * Tratamento seguro e fallback para **códigos de nucleotídeos degenerados IUPAC** (`R, Y, S, W, K, M, B, D, H, V, N`).
   * Exportação de relatórios JSON com encodificação **UTF-8 nativa (`ensure_ascii=False`)** e metadata de proveniência ISO 8601.
   * Suíte de testes unitários automatizados cobrindo 100% das funções centrais.

---

## 🔬 Arquitetura do Sistema e Módulos

```text
[Primers "Naturais" / Output NextLAMP]
                   │
                   ▼
┌──────────────────────────────────────────────────────────────┐
│  MÓDULO 1: Engenharia de Quimerização (chimerization.py)      │
│  - Suporte a ligantes nucleosídicos e não-nucleosídicos (HEG) │
│  - Matriz de Alinhamento Termodinâmico (SantaLucia 1998)      │
│  - Parser de chaves FIP/BIP fundidos e códigos IUPAC          │
│  - Superação do Art. 10, IX da LPI                            │
└──────────────────────────┬───────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  MÓDULO 2: Prova de Efeito Técnico (inventive_step.py)        │
│  - Ensaio Pareado: Primer Modificado vs. Controle Natural     │
│  - Fator de Supressão Termodinâmica por Equilíbrio Boltzmann  │
│  - Justificativa textual condicional sem contradições         │
│  - Demonstração de Atividade Inventiva (Art. 13 da LPI)       │
└──────────────────────────┬───────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  MÓDULO 3: Gerador do Pacote Legal (wipo_st26.py)             │
│  - XML WIPO ST.26 v1.3 com DTD estrito e feições misc_feature │
│  - SEQ ID NOs próprios para primers e ligantes isolados       │
│  - Minuta INPI Res. 248/2019 com Quadro Reivindicatório      │
└──────────────────────────┬───────────────────────────────────┘
```

---

## 📂 Estrutura de Arquivos

```text
patentlamp/
├── patentlamp/
│   ├── __init__.py            # Inicialização do pacote Python
│   ├── chimerization.py       # Módulo 1: Quimerização, Ligantes Sintéticos & Termodinâmica
│   ├── inventive_step.py      # Módulo 2: Prova de Efeito Técnico & Equilíbrio de Boltzmann
│   ├── wipo_st26.py           # Módulo 3: Gerador XML WIPO ST.26 & Minuta de Patente INPI
│   └── cli.py                 # Interface de Linha de Comando (CLI Multi-Set)
├── examples/                  # Casos de uso de exemplo e saídas validadas
├── tests/                     # Suíte de testes unitários automatizados (unittest)
├── setup.py                   # Script de instalação do pacote Python
├── requirements.txt           # Dependências do projeto
└── README.md                  # Este documento
```

---

## 💻 Guia de Uso

### 1. Execução via Linha de Comando (CLI)

Navegue até a raiz do repositório e invoque o **PatentLAMP** indicando o arquivo JSON gerado pelo **NextLAMP**:

```bash
export PYTHONPATH="${PWD}/patentlamp:${PYTHONPATH}"

python3 patentlamp/patentlamp/cli.py \
    --input-json results/babesia_paper_export/raw_data/nextlamp_babesia_results.json \
    --out-dir results/babesia_patent_package \
    --target-species "Babesia canis" \
    --applicant "Universidade Federal / LaBiOmicS" \
    --linker "HEG"
```

### 2. Parâmetros Disponíveis

| Parâmetro | Obrigatório | Descrição |
| :--- | :---: | :--- |
| `--input-json` | **Sim** | Caminho para o JSON contendo os primers desenhados pelo NextLAMP. |
| `--out-dir` | **Sim** | Diretório onde o pacote de patente completo será salvo. |
| `--target-species` | Não | Nome do organismo alvo (ex: `"Babesia canis"`, `"Toxoplasma gondii"`). Default: `"Babesia canis"`. |
| `--applicant` | Não | Nome do depositante/instituição para a listagem WIPO ST.26. |
| `--linker` | Não | Sequência ou ligante sintético para junção FIP/BIP (insensível a caixa). Opções padrão: `TAAA`, `TTTT`, `TTTTAAAATTTT`, `HEG`, `Spacer18`, `Spacer_C3`, `LNA`. Sequências customizadas (ex: `ACGTGA`) também são aceitas. Default: `"TAAA"`. |

---

## 📊 Arquivos Gerados no Pacote de Patente

Ao finalizar a execução, o **PatentLAMP** entrega os seguintes ativos na pasta de saída:

1. **`WIPO_ST26_Sequence_Listing.xml`**:
   * Arquivo XML rigorosamente estruturado no padrão WIPO ST.26 v1.3 com tags `INSDSeq_feature-table`, qualificadores `misc_feature` e entradas próprias de `SEQ ID NO` para primers e ligantes isolados, pronto para upload no e-Patentes (INPI), USPTO e EPO.
2. **`INPI_Patent_Specification_Draft.md`**:
   * Minuta do Relatório Descritivo em português contendo a contextualização da invenção, tabela de desempenho comparativo e o **Quadro Reivindicatório** completo (cautelado conforme a Resolução INPI nº 248/2019):
     * *Reivindicação 1 (Produto):* Conjunto de moléculas quiméricas sintéticas citando as SEQ ID NOs dos primers FIP/BIP e do ligante.
     * *Reivindicação 2 (Produto Complementar):* Primers externos F3, B3 e primers de alça LoopF/LoopB de forma dinâmica.
     * *Reivindicação 3 (Kit):* Composição do Kit Diagnóstico Isotérmico *In Vitro*.
     * *Reivindicação 4 (Processo):* Método de Diagnóstico *In Vitro*.
3. **`patentlamp_inventive_step_proofs.json`**:
   * Relatório estruturado em JSON com codificação UTF-8 nativa, metadata de proveniência ISO 8601, dados físico-químicos e métricas de $\Delta \Delta G$ e Fator de Supressão de Boltzmann.

---

## ⚖️ Conformidade Legal e Regulatória

* **Brasil (INPI):** Lei nº 9.279/1996 (LPI) - Arts. 10 (exclusões), 11 (novidade), 13 (atividade inventiva) e Resolução INPI nº 248/2019 (Biotecnologia).
* **Internacional (WIPO / PCT):** Padrão WIPO ST.26 para listagem de sequências em formato XML (obrigatório em todos os escritórios signatários da Convenção da União de Paris desde 01/07/2022).

---

## 👨‍💻 Autor & Licença

Desenvolvido pela equipe **LaBiOmicS / NextLAMP Project**.  
Licença de Software: MIT License / Academic & Commercial License.
