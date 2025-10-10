import sys
import os
from estrutura_proteica import Estrutura

def main():
    print("TESTE")

    estrutura = Estrutura(sys.argv[1])
    estrutura.resumo_estrutural()

if __name__ == "__main__":
    main()