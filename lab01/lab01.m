% Изучение дискретизации сигналов
function lab01
    n = input('Введите количество точек: ');
    step = input('Введите шаг: ');
    
    xMax = step * (n - 1) / 2;
    xArr = -xMax:step:xMax;
    x0Arr = -xMax:0.01:xMax;

    % Исходный сигнал
    yRec0Arr = rectangularPulse(x0Arr);
    yGau0Arr = gaussSignal(x0Arr);

    % Дискретизация
    yRecArr = rectangularPulse(xArr);
    yGauArr = gaussSignal(xArr);

    % Восстановление сигнала по формуле Котельникова
    yRecRArr = zeros(1, length(x0Arr));
    yGauRArr = zeros(1, length(x0Arr));

    for i = 1:length(x0Arr)
        for j = 1:n
            temp = normSinc((x0Arr(i) - xArr(j)) / step);
            yRecRArr(i) = yRecRArr(i) + yRecArr(j) * temp;
            yGauRArr(i) = yGauRArr(i) + yGauArr(j) * temp;
        end
    end

    figure;

    subplot(2,1,1);
    title('Прямоугольный импульс');
    hold on;
    grid on;
    plot(x0Arr, yRec0Arr, 'b');
    plot(x0Arr, yRecRArr, 'k');
    plot(xArr, yRecArr, '.r');
    legend('Исходный', 'Восстановленный', 'Дискретный');

    subplot(2,1,2);
    title('Сигнал Гаусса');
    hold on;
    grid on;
    plot(x0Arr, yGau0Arr, 'b');
    plot(x0Arr, yGauRArr, 'k');
    plot(xArr, yGauArr, '.r');
    legend('Исходный', 'Восстановленный', 'Дискретный');

end

function [y] = rectangularPulse(x)
    c = 3;
    y = zeros(size(x));

    y(abs(x) < c) = 1;
end

function [y] = gaussSignal(x)
    A = 1;
    sigma = 5;

    y = A * exp(-(x / sigma).^2);
end

function [y] = normSinc(x)
    if x ~= 0
        y = sin(pi * x) / (pi * x);
    else
        y = 1;
    end
end



    
    