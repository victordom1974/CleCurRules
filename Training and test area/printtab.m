function mess = printtab(tab, varargin)
    mess = "";  % Inicializar string vacío
    
    % Generar la tabla en formato LaTeX
    for i = 1:size(tab,1)
        % Primera columna en entero
        line = sprintf('%5d', tab(i,1));  
        
        % Resto de columnas en formato científico con 3 cifras significativas
        for j = 2:size(tab,2)
            line = strcat(line, sprintf(' & %.2e', tab(i,j)));
        end
        
        % Agregar doble barra al final de la fila
        line = strcat(line, ' \\\\\n');
        
        % Guardar en variable `mess`
        mess = strcat(mess, line);
    end

    % Si no se especifica un archivo de salida, mostrar en pantalla
    if nargin == 1 || isempty(varargin{1})
        sprintf(mess)
    else
        % Guardar en archivo especificado
        output_file = varargin{1};
        fid = fopen(output_file, 'w');  
        
        if fid == -1
            error('No se pudo abrir el archivo %s.', output_file);
        end
        
        fwrite(fid, mess);
        fclose(fid);
    end
end
