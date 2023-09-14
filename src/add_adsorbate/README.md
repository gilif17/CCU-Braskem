# Add Adsorbate to POSCAR Script

## Descrição

Este script Python permite adicionar um adsorvato a uma superfície em todas as posições possíveis de um tipo de site especificado. O script aceita cinco argumentos de linha de comando: o caminho para o arquivo POSCAR de entrada, o caminho para o arquivo POSCAR de saída, o símbolo do átomo adsorvato, o tipo de site e a altura acima do site para colocar o adsorvato.

## Requisitos

- Python 3.x
- ASE (Atomic Simulation Environment)
- SciPy
- NumPy

Você pode instalar os pacotes necessários usando o seguinte comando:

```bash
pip install -r requirements
```

## Uso

Para usar este script, você precisa especificar os seguintes argumentos de linha de comando:

- `input_poscar`: O caminho para o arquivo POSCAR de entrada que representa a superfície.
- `output_poscar`: O caminho para o arquivo POSCAR de saída onde a nova estrutura será escrita.
- `adsorbate_atom`: O símbolo do átomo adsorvato (por exemplo, 'H').
- `site`: O tipo de site para adicionar o adsorvato ('top', 'bridge' ou 'hole').
- `height`: A altura acima do site para colocar o adsorvato.

### Exemplo de Uso

```bash
python add_adsorbate.py input_poscar_path output_poscar_path H top 1.0
```

Neste comando:

- `input_poscar_path` é o caminho para o arquivo POSCAR de entrada.
- `output_poscar_path` é o caminho para o arquivo POSCAR de saída.
- `H` é o símbolo do átomo adsorvato.
- `top` é o tipo de site.
- `1.0` é a altura.