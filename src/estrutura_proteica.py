from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
from relatorios import Ramachandran, MapaContato, PontesHidrogenio, Sasa, Conformacao3d

class Estrutura:
    def __init__(self, arquivo):

        self.arquivo = arquivo
        try:
            self.estrutura, self.cabecalho = self.__carregamento_arquivo_proteina()
            self.modelos = list(self.estrutura.get_models())

            self.cadeias = []
            self.residuos = []
            self.atomos = []
            self.sequencias_polipeptidica = []
            self.__extracao_dados()
            
            sequencia_polipep_completa = ''.join(self.sequencias_polipeptidica)
            self.massa_molecular = self.__calculo_massa_molecular(sequencia_polipep_completa)
            self.ponto_isoeletrico = self.__calculo_ponto_isoletrico(sequencia_polipep_completa)

            destino_arquivos = "../results"
            self._gerar_graficos(destino_arquivos)

        except Exception as e:
            raise ValueError(f"* ERRO: Extração de metadados da proteína {e}.")

    def __carregamento_arquivo_proteina(self):
        try:
            parser = PDBParser(QUIET=True)
            estrutura = parser.get_structure('protein', self.arquivo)
            cabecalho = parser.get_header()
            return estrutura, cabecalho
        
        except FileNotFoundError:
            raise FileNotFoundError(f"ERRO: Arquvio da proteína não encontrado {self.arquivo}.")
        except Exception:
            raise ValueError(f"ERRO: Carregamento de arquivo da estrutura da proteína {self.arquivo}.")
        
    def __extracao_dados(self):
        pp_builder = PPBuilder()
        for modelo in self.modelos:
            for cadeia in modelo:
                self.cadeias.append(cadeia)
                polipeptideos = pp_builder.build_peptides(cadeia)
                for polipeptideo in polipeptideos:
                    self.sequencias_polipeptidica.append(str(polipeptideo.get_sequence()))
                for residuo in cadeia:
                    self.residuos.append(residuo)
                    for atomo in residuo:
                        self.atomos.append(atomo)
        
    def __calculo_massa_molecular(self, sequencia_polipep_completa):
        try:
            massa_molecular = molecular_weight(sequencia_polipep_completa, "protein")

        except Exception:
            massa_molecular = len(sequencia_polipep_completa) * 110
            print(f"ERRO: Cálculo da massa molecular. Será atribuido pelo tamanho da sequencia * {110}")

        return massa_molecular / 1000
    
    def __calculo_ponto_isoletrico(self, sequencia_polipep_completa):
        try:
            iep = IP(sequencia_polipep_completa)
            pi = iep.pi()

        except Exception:
            pi = 6.5
            print(f"ERRO: Cálculo do ponto isoelétrico. Será atribuído o valor {pi}")
        
        return pi
    
    def _gerar_graficos(self, destino):
        self.ramachandran = Ramachandran(self, destino)
        self.mapa_contato = MapaContato(self, destino)
        self.pontes_H = PontesHidrogenio(self, destino)
        self.sasa = Sasa(self, destino)
        self.conformacao = Conformacao3d(self, destino)
