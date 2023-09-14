# Projeto CCU-Braskem

## Descrição

O projeto CCU-Braskem é uma iniciativa focada no desenvolvimento de ferramentas e soluções para simulações de dinâmica molecular e análises estruturais em nível atômico. Este repositório hospeda uma série de scripts e utilitários Python desenvolvidos como parte desta iniciativa.

## Estrutura do Repositório

O repositório está estruturado da seguinte maneira:

- `src/`: Contém todos os scripts Python desenvolvidos para o projeto.
- `data/`: Diretório para armazenar arquivos de dados usados e gerados pelos scripts.
- `docs/`: Contém documentação detalhada para cada script, incluindo descrições, requisitos e exemplos de uso.
- `tests/`: Diretório que contém testes unitários para validar a funcionalidade dos scripts.

## Requisitos

Para executar os scripts disponíveis neste repositório, você precisará ter Python 3.x instalado, junto com várias bibliotecas Python, incluindo ASE, SciPy e NumPy. Você pode instalar todas as dependências necessárias usando o seguinte comando:

```bash
pip install -r requirements.txt
```

## Uso

Cada script no diretório `src/` é um utilitário independente que pode ser executado a partir da linha de comando. Consulte a documentação individual de cada script ou os cabeçalhos em cada diretório para obter instruções detalhadas sobre como usar cada script.

### Exemplo de Uso

Para usar o script `add_adsorbate.py`, você usaria um comando como o seguinte:

```
cd src/add_adsorbate
python src/add_adsorbate.py input_poscar_path output_poscar_path H top 1.0
```

## Contribuindo

Estamos abertos a contribuições para este projeto. Se você deseja contribuir, por favor:

1. Faça um fork do repositório.
2. Crie uma nova branch para sua feature ou correção de bug.
3. Faça suas alterações e teste-as completamente.
4. Envie um pull request detalhando as alterações feitas.

## Contato

Para qualquer dúvida ou feedback, por favor, entre em contato com a equipe de desenvolvimento através do [email](mailto:joaocassianox7x@gmail.com).

## Licença

Este projeto está licenciado sob a licença MIT. Veja o arquivo `LICENSE` para mais detalhes.