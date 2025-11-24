# PDBAnalyzer

Software para análise e geração de relatórios de estruturas PDB introdutórios.

## Descrição

Esta aplicação realiza a leitura de dados provenientes de um arquivo em formato PDB. Por meio do mesmo, cria gráficos de: Ramachandran, ligações de hidrogênio,
mapa de contatos, distruição SASA e de disposição tridimensional de carbonos alfa. Juntamente, gera um PDF integralizando estes artefatos e informações relacionadas.

## Instalação

### Pré requisitos

- Python 3.8
- pip

### Instalação das dependências

```bash
pip install -r requirements.txt
```

## Exemplo

Com o arquivo desejado e baixado, realize o comando abaixo com o primeiro argumento igual o exemplificado e, se o arquivo PDB estiver em /data, somente altere o restante
para o nome do arquivo desejado, caso contrário coloque o destino correto.

```bash
python src/pdb_analyzer.py data/1MDM.pdb
```

## Resultados

Os resultados das execuções serão destinados na pasta /results deste projeto, contendo cada gráfico em formato PNG e o PDF contendo todas as informações e imagens.

## Agradecimentos

- Programa Interunidades de Pós-Graduação em Bioinformática da UFMG - Universidade Federal de Minas Gerais. 