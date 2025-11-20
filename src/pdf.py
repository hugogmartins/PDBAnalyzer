import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import datetime

class GeradorPDF:
    def __init__(self, estrutura, destino):
        self.estrutura = estrutura
        self.destino = destino
        self.nome_molecula = self.estrutura.cabecalho['idcode']
        self.__gerar_pdf()

    def __gerar_pdf(self):
        try:
            caminho_pdf = os.path.join(self.destino, f"{self.nome_molecula}_relatorio_estrutural.pdf")

            with PdfPages(caminho_pdf) as pdf:
                self.__criar_capa(pdf)

                self.__adicionar_grafico(pdf, f'{self.nome_molecula}_ramachandran.png', "Gráfico de Ramachandran", "Distribuição dos ângulos dihedrais Phi e Psi dos resíduos protéicos.")
                self.__adicionar_grafico(pdf, f'{self.nome_molecula}_mapa_contato.png', "Mapa de Contatos", "Matriz de contatos entre carbonos alfa.")
                self.__adicionar_grafico(pdf, f'{self.nome_molecula}_ligacoes_hidrogenio.png', "Rede de Ligações de Hidrogênio", "Matriz e distribuição de ligações de hidrogênio na estrutura.")
                self.__adicionar_grafico(pdf, f'{self.nome_molecula}_distribuicao_sasa.png', "Análise de Acessibilidade", "Distribuição da Área Superficial Acessível ao Solvente.")
                self.__adicionar_grafico(pdf, f'{self.nome_molecula}_estrutura.png', "Gráfico 3D de posições de Carbonos Alfa", "")

                return caminho_pdf
        except Exception:
            print("ERRO: Geração de PDF.")
    
    def __criar_capa(self, pdf):
        fig, ax = plt.subplots(figsize=(8.27, 11.69))
        fig.patch.set_facecolor('white')
        ax.axis('off')

        ax.text(0.5, 0.9, f"Relatório PDBAnalyzer {self.nome_molecula}", ha='center', va='center', fontsize=16, fontweight='bold', transform=ax.transAxes)

        ax.text(0.1, 0.75, "Análise do arquivo:", ha='left', va='center', fontsize=12, fontweight='bold', transform=ax.transAxes)

        resumo_pdb = {
            'ID' : self.estrutura.cabecalho['idcode'],
            'Nome' : self.estrutura.cabecalho['name'],
            'Método' : self.estrutura.cabecalho['structure_method'],
            'Resolução' : str(self.estrutura.cabecalho['resolution']) + 'A',
            'Autor' : self.estrutura.cabecalho['author'],
            'Data de disponibilidade' : self.estrutura.cabecalho['release_date'],
            'Link' : "https://www.rcsb.org/structure/" + self.estrutura.cabecalho['idcode']
        }

        self.__criar_resumo_dados(ax, 0.7, resumo_pdb)

        ax.text(0.1, 0.45, "Análise estrutural:", ha='left', va='center', fontsize=12, fontweight='bold', transform=ax.transAxes)

        resumo_estrutural = {
            'Modelos': str(len(self.estrutura.modelos)) + (' unidade.' if len(self.estrutura.modelos) == 1 else ' unidades.'),
            'Cadeias': str(len(self.estrutura.cadeias)) + (' unidade.' if len(self.estrutura.cadeias) == 1 else ' unidades.'),
            'Resíduos': str(len(self.estrutura.residuos)) + (' unidade.' if len(self.estrutura.residuos) == 1 else ' unidades.'),
            'Átomos' : str(len(self.estrutura.atomos)) + (' unidade.' if len(self.estrutura.atomos) == 1 else ' unidades.'),
            'Massa molecular': round(self.estrutura.massa_molecular, 2),
            'Pontos isoelétricos': round(self.estrutura.ponto_isoeletrico, 2)
        }

        self.__criar_resumo_dados(ax, 0.4, resumo_estrutural)

        data = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        ax.text(0.5, 0.05, f'PDBAnalyzer - Data de geração: {data}', ha='center', va='center', fontsize=10, transform=ax.transAxes)

        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def __adicionar_grafico(self, pdf, nome_arquivo, titulo, descricao):
        caminho = os.path.join(self.destino, nome_arquivo)
        if os.path.exists(caminho):
            fig = plt.figure(figsize=(8.27, 11.69))
            fig.patch.set_facecolor('white')

            gs = plt.GridSpec(3, 1, height_ratios=[0.1, 0.7, 0.2])

            ax_titulo = plt.subplot(gs[0])
            ax_titulo.axis('off')
            ax_titulo.text(0.5, 0.5, titulo, ha='center', va='center', fontsize=14, fontweight='bold', transform=ax_titulo.transAxes)

            ax_grafico = plt.subplot(gs[1])
            img = plt.imread(caminho)
            ax_grafico.imshow(img)
            ax_grafico.axis('off')

            ax_desc = plt.subplot(gs[2])
            ax_desc.axis('off')
            ax_desc.text(0.5, 0.5, descricao, ha='center', va='center', fontsize=10, style='italic', wrap=True, transform=ax_desc.transAxes)

            data = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
            ax_desc.text(0.5, 0.3, f'PDBAnalyzer - Data de geração: {data}', ha='center', va='center', fontsize=10, transform=ax_desc.transAxes)

            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
        
        else:
            fig, ax = plt.subplots(figsize=(8.27, 11.69))
            fig.patch.set_facecolor('white')
            ax.axis('off')
            
            ax.text(0.5, 0.7, titulo, ha='center', va='center', fontsize=14, fontweight='bold', transform=ax.transAxes)
            
            ax.text(0.5, 0.5, f'Gráfico {nome_arquivo} não disponível', ha='center', va='center', fontsize=12, color='red', transform=ax.transAxes)
            
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

    def __criar_resumo_dados(self, eixo, eixo_y, dicionario_texto):
        if dicionario_texto:
            for chave, valor in dicionario_texto.items():
                texto = f"{chave}: {valor}"
                eixo_y = self.__tratamento_quebra_de_linha(eixo, eixo_y, texto, 10)
        return eixo_y

    def __tratamento_quebra_de_linha(self, eixo, eixo_y, texto, fonte):
        if (len(texto) * fonte * 0.75) > 500:
            palavras = texto.split(' ')
            linha = ""
            for palavra in palavras:
                linha = linha + " " + palavra if len(linha) != 0 else palavra
                largura = len(linha) * fonte * 0.75
                
                if largura > 500:
                    eixo.text(0.1, eixo_y, linha, ha='left', va='center', fontsize=fonte, transform=eixo.transAxes)
                    eixo_y = eixo_y - 0.025
                    linha = ""
            
            if len(linha) != 0:
                eixo.text(0.1, eixo_y, linha, ha='left', va='center', fontsize=fonte, transform=eixo.transAxes)

        else:
            eixo.text(0.1, eixo_y, texto, ha='left', va='center', fontsize=fonte, transform=eixo.transAxes)
        
        return eixo_y - 0.025