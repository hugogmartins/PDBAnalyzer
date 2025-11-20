from abc import ABC, abstractmethod
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.vectors import calc_dihedral
from Bio.PDB import ShrakeRupley
import numpy as np
import matplotlib.pyplot as plt
import os

class Relatorio(ABC):

    def __init__(self, estrutura, destino_arquivos):
        self.estrutura = estrutura
        self.destino_arquivos = destino_arquivos
        self.__integralizacao(destino_arquivos)

    @abstractmethod
    def _calcular(self) -> dict:
        pass

    @abstractmethod
    def _plot(self) -> str:
        pass

    def __integralizacao(self, destino_arquivos):
        self.resultado = self._calcular()
        destino_arquivos = self._plot(self.resultado, destino_arquivos)

class Ramachandran(Relatorio):

    def _calcular(self):
        try:
            angulos_phi = []
            angulos_psi = []
            residuos = self.estrutura.residuos
            for i in range(1, len(residuos) - 1):
                residuo_anterior = residuos[i - 1]
                residuo_atual = residuos[i]
                residuo_posterior = residuos[i + 1]

                if not (is_aa(residuo_anterior) and is_aa(residuo_atual) and is_aa(residuo_posterior)):
                    continue
                
                try:
                    phi = calc_dihedral(
                        residuo_anterior['C'].get_vector(),
                        residuo_atual['N'].get_vector(),
                        residuo_atual['CA'].get_vector(),
                        residuo_atual['C'].get_vector()
                    )
                    psi = calc_dihedral(
                        residuo_atual['N'].get_vector(),
                        residuo_atual['CA'].get_vector(),
                        residuo_atual['C'].get_vector(),
                        residuo_posterior['N'].get_vector()
                    )
                    phi_graus = np.degrees(phi)
                    psi_graus = np.degrees(psi)

                    angulos_phi.append(phi_graus)
                    angulos_psi.append(psi_graus)
                
                except KeyError:
                    continue

            resultado = {
                'phi': angulos_phi,
                'psi': angulos_psi
            }

            return resultado

        except Exception:
            print("ERRO: Cálculo dos ângulos de torção(Phi e Psi) para Ramachandran.")

    def _plot(self, resultado, destino_arquivo="ramachandran.png"):
        angulos_phi = resultado['phi']
        angulos_psi = resultado['psi']
        try:
            figura, eixo = plt.subplots(figsize=(10, 8))

            plt.scatter(angulos_phi, angulos_psi, alpha=0.6, s=20, c='blue', edgecolors='black', linewidths=0.5)
            
            plt.xlim(-180, 180)
            plt.ylim(-180, 180)

            plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
            plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
            plt.grid(True, alpha=0.3)

            plt.title("Ramachandran", fontsize=16, fontweight='bold')
            plt.xlabel("Ângulos PHI (°).", fontsize=12)
            plt.ylabel("Ângulos PSI (°)", fontsize=12)

            plt.text(-150, -150, "β-sheet", fontsize=12, color='red', alpha=0.7)
            plt.text(-100, 50, "α-helix", fontsize=12, color='red', alpha=0.7)
            plt.text(50, 50, "α-helix canhota", fontsize=12, color='red', alpha=0.7)

            os.makedirs(destino_arquivo, exist_ok=True)
            destino_arquivo = os.path.join(destino_arquivo, f"{self.estrutura.cabecalho['idcode']}_ramachandran.png")
            plt.tight_layout()
            plt.savefig(destino_arquivo, dpi=300, bbox_inches='tight')
            plt.close(figura)

            print(f"* Gráfico de Ramachandran salvo em: {destino_arquivo}")
            return destino_arquivo
        
        except Exception:
            raise ValueError("* ERRO: Criação do gráfico de Ramachandran.")
        
        
class MapaContato(Relatorio):
    def _calcular(self, limite_distancia=8.0):
        try:
            atomos_c_alfa = []

            residuos = self.estrutura.residuos

            for residuo in residuos:
                if is_aa(residuo):
                    try:
                        c_alfa = residuo['CA']
                        atomos_c_alfa.append(c_alfa.get_coord())
                    except KeyError:
                        continue
            
            total_residuos = len(atomos_c_alfa)

            if total_residuos == 0:
                raise ValueError("* ERRO: Nenhum Carbono Alfa encontrado na estrutura.")

            matriz_contato = np.zeros((total_residuos, total_residuos))

            for i in range(total_residuos):
                for j in range(i + 1, total_residuos):
                    distancia = np.linalg.norm(atomos_c_alfa[i] - atomos_c_alfa[j])
                    if distancia <= limite_distancia:
                        matriz_contato[i, j] = 1
                        matriz_contato[j, i] = 1
            
            resultado = {
                'matriz': matriz_contato,
                'carbonos_alfa': atomos_c_alfa
            }

            return resultado
        
        except Exception:
            raise ValueError("* ERRO: Cálculo do mapa de contatos.")
        
    def _plot(self, resultado, destino_arquivo="mapa_contato.png"):
        matriz = resultado['matriz']
        tamanho_matriz = matriz.shape[0]

        figura, eixos = plt.subplots(figsize=(12,10))
        imagem = eixos.imshow(matriz, cmap='binary', interpolation='nearest', origin='lower', vmin=0, vmax=1)

        barrac = plt.colorbar(imagem, ax=eixos)
        barrac.set_label("Contato: 1 - Presente ; 2 - Ausente")
        barrac.set_ticks([0, 1])
        barrac.set_ticklabels(['0', '1'])

        eixos.set_title("Mapa de Contatos", fontsize=16, fontweight='bold')
        eixos.set_xlabel("Número do Resíduo.", fontsize=12)
        eixos.set_ylabel("Número do Resíduo.", fontsize=12)
        
        if tamanho_matriz > 0:
            if tamanho_matriz <= 50:
                eixos.set_xticks(range(tamanho_matriz))
                eixos.set_yticks(range(tamanho_matriz))
                eixos.set_xticklabels(range(tamanho_matriz), rotation=90, fontsize=8)
                eixos.set_yticklabels(range(tamanho_matriz), fontsize=8)
            
            else:
                passo = max(1, tamanho_matriz // 20)
                indices = range(0, tamanho_matriz, passo)
                eixos.set_xticks(indices)
                eixos.set_yticks(indices)
                eixos.set_xticklabels(indices, rotation=90, fontsize=8)
                eixos.set_yticklabels(indices, fontsize=8)
            
        eixos.grid(True, which='both', color='lightgray', linestyle='-', linewidth=0.5, alpha=0.3)

        os.makedirs(destino_arquivo, exist_ok=True)
        destino_arquivo = os.path.join(destino_arquivo, f"{self.estrutura.cabecalho['idcode']}_mapa_contato.png")
        plt.tight_layout()
        plt.savefig(destino_arquivo, dpi=300, bbox_inches='tight')
        plt.close(figura)

        print(f"* Mapa de Contatos salvo em: {destino_arquivo}")
        return destino_arquivo
    
class PontesHidrogenio(Relatorio):

    def _calcular(self):
        resultado = self.__encontrar_pontes()
        resultado = self.__criar_matriz(resultado)
        return resultado
    
    def __encontrar_pontes(self, limite_distancia=3.5):
        try:
            pontes = []
            doadores = ['N', 'OG', 'OG1', 'NE', 'NH1', 'NH2', 'ND2', 'NH2', 'NZ']
            receptores = ['O', 'OD1', 'OD2', 'OE1', 'OE2', 'OH']
            cadeias = self.estrutura.cadeias

            atomos = []
            informacoes_residuos = {}

            for cadeia in cadeias:
                for residuo in cadeia:
                    if is_aa(residuo):
                        residuo_id = f"{cadeia.id}:{residuo.resname}{residuo.id[1]}"
                        informacoes_residuos[residuo_id] = {
                            'cadeia': cadeia.id,
                            'residuo': f"{residuo.resname} {residuo.id[1]}",
                            'residuo_numero': residuo.id[1],
                            'residuo_nome': residuo.resname
                        }

                        for atomo_nome in doadores:
                            try:
                                atomo = residuo[atomo_nome]
                                atomos.append({
                                    'atomo': atomo,
                                    'residuo_id': residuo_id,
                                    'tipo': 'doador',
                                    'atomo_nome': atomo_nome 
                                })
                            except KeyError:
                                continue

                        for atomo_nome in receptores:
                            try:
                                atomo = residuo[atomo_nome]
                                atomos.append({
                                    'atomo': atomo,
                                    'residuo_id': residuo_id,
                                    'tipo': 'receptor',
                                    'atomo_nome': atomo_nome 
                                })
                            except KeyError:
                                continue
            
            for i, atomo_i in enumerate(atomos):
                for j, atomo_j in enumerate(atomos[i + 1:], i + 1):
                    if atomo_i['tipo'] == atomo_j['tipo']:
                        continue
                    distancia = atomo_i['atomo'] - atomo_j['atomo']

                    if distancia <= limite_distancia:
                        if atomo_i['tipo'] == 'doador':
                            doador = atomo_i
                            receptor = atomo_j
                        else:
                            doador = atomo_j
                            receptor = atomo_i
                        
                        pontes.append({
                            'residuo_doador': doador['residuo_id'],
                            'residuo_receptor': receptor['residuo_id'],
                            'atomo_doador': doador['atomo_nome'],
                            'atomo_receptor': receptor['atomo_nome'],
                            'distancia': distancia,
                            'cadeia_doador': informacoes_residuos[doador['residuo_id']]['cadeia'],
                            'cadeia_receptor': informacoes_residuos[receptor['residuo_id']]['cadeia'],
                            'doador_id': informacoes_residuos[doador['residuo_id']]['residuo_numero'],
                            'receptor_id': informacoes_residuos[receptor['residuo_id']]['residuo_numero']
                        })
            
            resultado = {
                'pontes': pontes,
                'informacoes': informacoes_residuos
            }

            return resultado
        
        except Exception:
            raise ValueError("* ERRO: Identificação de pontes de Hidrogênio.")
    
    def __criar_matriz(self, resultado):
        try:
            pontes = resultado['pontes']

            residuos = sorted(set([ponte['residuo_doador'] for ponte in pontes] + [ponte['residuo_receptor'] for ponte in pontes]))
            residuo_idx = {res: idx for idx, res in enumerate(residuos)}
            total_residuos = len(residuos)

            pontes_matriz = np.zeros((total_residuos, total_residuos))

            for ponte in pontes:
                i = residuo_idx[ponte['residuo_doador']]
                j = residuo_idx[ponte['residuo_receptor']]
                pontes_matriz[i, j] = 1
                pontes_matriz[j, i] = 1

            informacoes_residuos_ordenadas = []
            for residuo in residuos:
                partes = residuo.split(':')
                if len(partes) == 2:
                    cadeia, rest = partes
                    residuo_nome = rest[:3]
                    residuo_num_str = rest[3:]
                    try:
                        residuo_numero = int(residuo_num_str)
                        informacoes_residuos_ordenadas.append({
                            'cadeia': cadeia,
                            'residuo': f"{residuo_nome} {residuo_numero}",
                            'residuo_numero': residuo_numero,
                            'residuo_nome': residuo_nome
                        })
                    except ValueError:
                        continue
            
            resultado ={
                'matriz': pontes_matriz,
                'informacoes_ordenadas': informacoes_residuos_ordenadas
            }
            
            return resultado
        
        except Exception as e:
            raise ValueError(f"Erro ao criar matriz de ligações de H.")

    def _plot(self, resultado, destino_arquivo="ligacoes_hidrogenio.png"):
        try:
            informacoes_residuos = resultado['informacoes_ordenadas']
            matriz = resultado['matriz']
            figura, (eixo1, eixo2) = plt.subplots(1, 2, figsize=(16, 8))

            total_residues = len(informacoes_residuos)

            im = eixo1.imshow(matriz, cmap='Reds', interpolation='nearest', origin='lower')
            plt.colorbar(im, ax=eixo1, label="Ponte de Hidrogênio (1 = Presente, 0 = Ausente)")

            eixo1.set_title("Matriz de ligações de Hidrogênio", fontsize=14, fontweight='bold')
            eixo1.set_xlabel("Índice do Resíduo")
            eixo1.set_ylabel("Índice do Resíduo")
            eixo1.grid(True, alpha=0.3)

            if total_residues <= 50:
                rotulos_residuos = [f"{info['residuo_nome']}{info['residuo_numero']}" for info in informacoes_residuos]
                eixo1.set_xticks(range(total_residues))
                eixo1.set_yticks(range(total_residues))
                eixo1.set_xticklabels(rotulos_residuos, rotation=90, fontsize=6)
                eixo1.set_yticklabels(rotulos_residuos, fontsize=6)
            
            pontes_por_residuo = np.sum(matriz, axis=1)

            barras = eixo2.bar(range(total_residues), pontes_por_residuo, color='lightcoral', edgecolor='darkred', alpha=0.7)
            eixo2.set_xlabel("Índice do Resíduo")
            eixo2.set_ylabel("Número de ligações de Hidrogênio")
            eixo2.set_title("Distribuição de ligações de H por Resíduo", fontsize=14, fontweight='bold')
            eixo2.grid(True, alpha=0.2)

            if total_residues > 0:
                max_pontes = np.max(pontes_por_residuo)
                max_indices = np.where(pontes_por_residuo == max_pontes)[0]
                for idx in max_indices:
                    barras[idx].set_color('red')
            
            os.makedirs(destino_arquivo, exist_ok=True)
            destino_arquivo = os.path.join(destino_arquivo, f"{self.estrutura.cabecalho['idcode']}_ligacoes_hidrogenio.png")
            plt.tight_layout()
            plt.savefig(destino_arquivo, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"Rede de ligações de H salva como: {destino_arquivo}")
            return destino_arquivo
        
        except Exception as e:
            raise ValueError(f"Erro ao gerar rede de ligações de H.")
        
class Sasa(Relatorio):

    def _calcular(self):
        resultado = self.__obter_sasa()
        resultado = self.__classificar_acessibilidade(resultado)
        return resultado
        
    def __obter_sasa(self, circuferencia_sonda=1.4, n_pontos=100):
        try:
            shrake_rupley = ShrakeRupley(probe_radius=circuferencia_sonda, n_points=n_pontos)
            shrake_rupley.compute(self.estrutura.estrutura, level="R")

            dados_sasa = []
            total_sasa = 0

            cadeias = self.estrutura.cadeias

            for cadeia in cadeias:
                for residuo in cadeia:
                    if is_aa(residuo):
                        residuo_sasa = residuo.sasa
                        total_sasa += residuo_sasa

                        dados_sasa.append({
                            'cadeia': cadeia.id,
                            'residuo': f"{residuo.resname} {residuo.id[1]}",
                            'residuo_numero': residuo.id[1],
                            'residuo_nome': residuo.resname,
                            'sasa': residuo_sasa
                        })
            
            resultado = {
                'dados': dados_sasa,
                'total': total_sasa
            }
            return resultado
        
        except Exception:
            raise ValueError(f"Erro ao calcular SASA.")
    
    def __classificar_acessibilidade(self, resultado, metodo='percentil'):
        try:
            dados = resultado['dados']
            valores_sasa = [dado['sasa'] for dado in dados]
            if metodo == 'percentil':
                percentil33 = np.percentile(valores_sasa, 33)
                percentil66 = np.percentile(valores_sasa, 66)
                
                for dado in dados:
                    if dado['sasa'] <= percentil33:
                        dado['acessibilidade'] = 'Inacessível'
                        dado['classe'] = 0
                    elif dado['sasa'] <= percentil66:
                        dado['acessibilidade'] = 'Parcialmente Exposto'
                        dado['classe'] = 1
                    else:
                        dado['acessibilidade'] = 'Exposto'
                        dado['classe'] = 2
            return dados

        except Exception:
            raise ValueError
        
    def _plot(self, resultado, destino_arquivo="distribuicao_sasa.png"):

        try:
            figura, eixos = plt.subplots(2, 2, figsize=(15, 12))
            valores_sasa = [dado['sasa'] for dado in resultado]
            residuos_numeros = [dado['residuo_numero'] for dado in resultado]

            eixos[0, 0].plot(residuos_numeros, valores_sasa, 'b-', linewidth=1, alpha=0.7)
            eixos[0, 0].fill_between(residuos_numeros, valores_sasa, alpha=0.3)
            eixos[0, 0].set_title("Distribuição de valores SASA.", fontweight='bold')
            eixos[0, 0].set_xlabel("Número do Resíduo.")
            eixos[0, 0].set_ylabel("SASA Å²")
            eixos[0, 0].grid(True, alpha=0.3)

            eixos[0, 1].hist(valores_sasa, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
            eixos[0, 1].axvline(np.mean(valores_sasa), color='red', linestyle='--', label=f"Média {np.mean(valores_sasa):.2f} Å²")
            eixos[0, 1].axvline(np.median(valores_sasa), color='green', linestyle='--', label=f"Mediana {np.median(valores_sasa):.2f} Å²")
            eixos[0, 1].set_xlabel("SASA Å²")
            eixos[0, 1].set_ylabel("Frequência")
            eixos[0, 1].set_title("Distribuição de valores de SASA.", fontweight='bold')
            eixos[0, 1].legend()
            eixos[0, 1].grid(True, alpha=0.3)

            if 'acessibilidade' in resultado[0]:
                contador_acessibilidade = {}
                for data in resultado:
                    acessibilidade = data['acessibilidade']
                    contador_acessibilidade[acessibilidade] = contador_acessibilidade.get(acessibilidade, 0) + 1
                cores = ['lightcoral', 'gold', 'lightgreen']
                eixos[1, 0].pie(contador_acessibilidade.values(), labels=contador_acessibilidade.keys(), autopct='%1.1f%%', colors=cores, startangle=90)
                eixos[1, 0].set_title("Classificação de Acessibilidade dos Resíduos.", fontweight='bold')
            
            residuos_tipos = {}
            for data in resultado:
                residuo_nome = data['residuo_nome']
                if residuo_nome == data['residuo_nome']:
                    residuos_tipos[residuo_nome] = []
                residuos_tipos[residuo_nome].append(data['sasa'])

            residuo_media = {res: np.mean(sasas) for res, sasas in residuos_tipos.items()}
            residuos_organizados = sorted(residuo_media.items(), key=lambda x: x[1])

            residuos_nomes = [res[0] for res in residuos_organizados]
            residuos_medias = [res[1] for res in residuos_organizados]

            barras = eixos[1, 1].bar(range(len(residuos_nomes)), residuos_medias, color='lightblue', edgecolor='black')
            eixos[1, 1].set_xlabel("Tipo de resíduo")
            eixos[1, 1].set_ylabel("SASA Média Å²")
            eixos[1, 1].set_title("SASA Médio por tipo de Resíduo.", fontweight='bold')
            eixos[1, 1].set_xticks(range(len(residuos_nomes)))
            eixos[1, 1].set_xticklabels(residuos_nomes, rotation=45, ha='right')

            for barra, valor in zip(barras, residuos_medias):
                eixos[1, 1].text(barra.get_x() + barra.get_width() / 2, barra.get_height() + 0.5, f"{valor:.1f}", ha='center', fontsize=8)
            
            os.makedirs(destino_arquivo, exist_ok=True)
            destino_arquivo = os.path.join(destino_arquivo, f"{self.estrutura.cabecalho['idcode']}_distribuicao_sasa.png")
            plt.tight_layout()
            plt.savefig(destino_arquivo, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"Gráficos SASA salvos como: {destino_arquivo}")
            return destino_arquivo

        except Exception as e:
            raise ValueError(f"Erro ao gerar os gráficos SASA.")
    
class Conformacao3d(Relatorio):
    def _calcular(self):
        resultado = {}
        return resultado

    def _plot(self, resultado, destino):
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')

        for modelo in self.estrutura.estrutura:
            for cadeia in modelo:
                for residuo in cadeia:
                    try:
                        ca = residuo['CA']
                        x, y, z = ca.get_coord()
                        ax.scatter(x, y, z, c='blue', s=20)
                    except KeyError:
                        continue
        
        ax.set_title(f"Estrutura 3D dos carbonos alfa")
        ax.set_xlabel('x (A)')
        ax.set_ylabel('y (A)')
        ax.set_zlabel('z (A)')

        destino_arquivo = os.path.join(destino, f"{self.estrutura.cabecalho['idcode']}_estrutura.png")
        plt.savefig(destino_arquivo, dpi=150, bbox_inches='tight')
        plt.close()

        return destino_arquivo
