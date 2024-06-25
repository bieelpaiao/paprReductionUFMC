function [xPadded, yPadded] = padMatrices(x, y)
    sizes = [size(x) size(y)];
    higherSize = max(sizes);

    powers2 = 1:15;
    squareSizes = 2.^powers2;

    for i=1:length(squareSizes)
        if squareSizes(1,i) > higherSize
            linhas_desejadas = squareSizes(1,i);
            colunas_desejadas = squareSizes(1,i);
            break
        end
    end

    % Calcula o número de linhas e colunas a serem adicionadas
    linhas_adicionais_x = linhas_desejadas - size(x, 1);
    colunas_adicionais_x = colunas_desejadas - size(x, 2);

    linhas_adicionais_y = linhas_desejadas - size(y, 1);
    colunas_adicionais_y = colunas_desejadas - size(y, 2);
    
    % Preenche a matriz original com zeros para alcançar as dimensões desejadas
    xPadded = padarray(x, [linhas_adicionais_x, colunas_adicionais_x], 0, 'post');
    yPadded = padarray(y, [linhas_adicionais_y, colunas_adicionais_y], 0, 'post');
end