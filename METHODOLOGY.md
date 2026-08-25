# 🔬 Metodologia Científica: Desenho de Primers LAMP de Escala Genômica com o NextLAMP

Este documento detalha a **metodologia experimental e bioinformática** para o desenho de primers de amplificação isotérmica mediada por loop (**LAMP**) e primers de **Loop** (LoopF e LoopB) direcionados ao gênero *_Babesia_* com filtragem de especificidade tripla contra hospedeiro, vetores e protozoários relacionados.

---

## 🎯 1. Escopo e Objetivos do Ensaio
O ensaio LAMP visa a detecção sensível e específica de espécies do gênero *_Babesia_* (*B. canis*, *B. vogeli*, *B. gibsoni*, etc.), permitindo o diagnóstico em tempo real tanto em amostras clínicas do hospedeiro vertebrado quanto em vigilância entomológica/epidemiológica direta nos carrapatos vetores.

---

## 🛡️ 2. Arquitetura do Banco de Dados de Seleção Negativa (Background)

Para garantir **zero reação cruzada** em ensaios *Point-of-Care* de campo ou extrações diretas de sangue e vetores, o banco de dados de exclusão (*background*) do NextLAMP foi estruturado em **3 pilares biológicos**:

### Pillar A: Protozoários Filogeneticamente Próximos (Apicomplexa não-*Babesia*)
Evita falsos-positivos por semelhança evolutiva com outros coccídeos e hemoparasitas:
- *Theileria* spp.
- *Plasmodium* spp.
- *Toxoplasma gondii*
- *Neospora caninum*
- *Cryptosporidium* spp.
- *Eimeria* spp.

### Pillar B: Hospedeiro Vertebrado
Evita ligação não-específica com a vasta quantidade de DNA genômico do hospedeiro presente nas amostras clínicas de sangue:
- *Canis lupus familiaris* (Cão doméstico - 20 montagens de genoma do NCBI).

### Pillar C: Carrapatos Vetores (Arthropoda: Ixodidae)
Evita reação cruzada quando o teste for aplicado diretamente em homogeneizados de carrapatos coletados no campo ou na pele dos animais:
- ***Rhipicephalus sanguineus*** (Carrapato marrom do cão - TaxID: `7469`)
- ***Dermacentor reticulatus*** (Vetor de *B. canis* na Europa - TaxID: `160490`)
- ***Rhipicephalus microplus*** (Carrapato do boi - TaxID: `6941`)
- ***Ixodes scapularis*** (Carrapato ixodídeo - TaxID: `6945`)

---

## ⚙️ 3. Parâmetros Termodinâmicos e Distâncias Espaciais

O algoritmo **NextLAMP** avalia simultaneamente 8 oligonucleotídeos por conjunto (F3, F2, F1c, LoopF, B1c, LoopB, B2, B3) sob as seguintes restrições:

| Parâmetro / Oligonucletídeo | Faixa Aceita / Limiar | Função Biológica |
| :--- | :---: | :--- |
| **Temperatura de Anelamento (\(Tm\))** | \(55.0^\circ\text{C} - 68.0^\circ\text{C}\) | Faixa ideal para reação isotérmica com *Bst* Polymerase |
| **Conteúdo de GC (\(GC\%\))** | \(30.0\% - 70.0\%\) | Estabilidade da dupla fita |
| **Balanço Termodinâmico (\(Tm_{\text{balance}}\))** | Minimizar \(|Tm(F2)-Tm(B2)| + |Tm(F3)-Tm(B3)|\) | Amplificação simétrica e eficiente das duas fitas |
| **Distância F3 → F2** | \(0 - 20\text{ bp}\) | Região de clivagem e deslocamento da fita externa |
| **Distância F2 → F1c** | \(40 - 60\text{ bp}\) | Formação da alça da estrutura *dumbbell* |
| **Tamanho do Amplicon Interno (F2 → B2)** | \(120 - 180\text{ bp}\) | Região central de amplificação isotérmica |
| **Distância B1c → B2** | \(40 - 60\text{ bp}\) | Formação da alça direita do *dumbbell* |
| **Distância B2 → B3** | \(0 - 20\text{ bp}\) | Deslocamento da fita externa direita |
| **Primers de Loop (LoopF / LoopB)** | Regiões entre F1-F2 e B1-B2 | Aceleração da reação (amplificação em < 20 min) |

---

## ⚡ 4. Mecanismo de Alinhamento e Proveniência FAIR

1. **Geração de Candidatos:** Varredura em janela deslizante no genoma modelo de *_Babesia_* identificando milhões de sequências viáveis.
2. **Filtragem de Especificidade via Bowtie 2:** Alinhamento paralelo contra o banco completo de 16.2 GB (Apicomplexa + Cão + Carrapatos) descartando candidatos com pareamentos imprecisos ou não-específicos.
3. **Deduplicação de Locus:** Eliminação de conjuntos redundantes para garantir cobertura espacial ampla de marcadores genômicos únicos.
4. **Relatório FAIR:** Exportação em JSON (hash SHA-256 do genoma alvo e parâmetros), TSV (pronto para pedido em laboratório) e TXT (interpretação rápida).
