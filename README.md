# TCC - Comparação de performance de técnicas de DoA
Nome: Rafael Pintar Alevato

Orientador: Richard Demo Souza

Coorientador: Pedro Lemos

Instituição: UFSC (Universidade Federal de Satana Catarina)

# Dependências

* Biblioteca Eigen: `git submodule update --init --recursive`

* Simplicity Commander: `./scripts/install_commander`

# Para rodar a aplicação

Após a instalação do Commander:

1. Conecte apenas o kit da âncora em uma porta USB e execute
    ```
    ./scripts/flash.sh -f anchor_firmware/bin/anchor-firmware.s37
    ./scripts/flash.sh -f anchor_firmware/bin/anchor-bootloader.s37
    ```

2. Desconecte o kit da âncora e conecte apenas o kit da tag em uma porta USB e execute
    ```
    ./scripts/flash.sh -f tag_firmware/bin/tag-firmware.s37
    ```

3. Desconecte o kit da tag e conecte novamente o Kit da âncora e verifique qual porta VCOM foi conectado (Ex.: /dev/ttyACM0).

4. Compile e execute a aplicação principal
    ```
    cd main_app/
    make
    cd exe
    ./locator-host /dev/ttyACM0
    ```

    Se o firmware da âncora tiver sido gravado corretamente você deve ver a saída:
    ```
    AoX library init...
    BGLIB init...
    uart_port: /dev/ttyACM0
    Starting up...
    Resetting NCP target...
    System booted. Looking for tags...
    ```

5. Conecte novamente a tag ou alimente com bateria. Se tudo estiver correto você deve ver a saída dos ângulos processados
    ```
    azimuth: 130.970825, elevation: 132.902817
    azimuth: 137.350067, elevation: 184.111053
    azimuth: 122.667847, elevation: 130.076935
    azimuth: 116.089989, elevation: 117.508133
    azimuth: 119.015411, elevation: 125.327667
    ```
