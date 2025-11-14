import sys
from estrutura_proteica import Estrutura
from pdf import GeradorPDF

def main():

    estrutura = Estrutura(sys.argv[1])
    pdf = GeradorPDF(estrutura, "results")

if __name__ == "__main__":
    main()