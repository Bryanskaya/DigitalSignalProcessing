function lab02
    n = input('Введите количество точек: ');
    step = input('Введите шаг: ');
    
    xMax = step * (n - 1) / 2;
    
    xArr = -xMax:step:xMax;
    xNum = 0:length(xArr) - 1;

    % Исходный сигнал
    yRecArr = rectangularPulse(xArr);
    yGauArr = gaussSignal(xArr);

    % БПФ (с эффектом близнецов)
    tic
    yRecFFT = fft(yRecArr);
    yGauFFT = fft(yGauArr);
    toc

    % Устранение эффекта близнецов
    yRecFFT_NT = fftshift(yRecFFT);
    yGauFFT_NT = fftshift(yGauFFT);

    % ДПФ (с эффектом близнецов)
    tic
    yRecDFT = dft(yRecArr);
    yGauDFT = dft(yGauArr);
    toc

    % Устранение эффекта близнецов
    yRecDFT_NT = fftshift(yRecDFT);
    yGauDFT_NT = fftshift(yGauDFT);

    figure(1);
    drawFigure(1, xNum, yRecFFT, yRecFFT_NT, 'БПФ', 'прямоугольный сигнал');
    drawFigure(2, xNum, yGauFFT, yGauFFT_NT, 'БПФ', 'Гаусс');
    drawFigure(3, xNum, yRecDFT, yRecDFT_NT, 'ДПФ', 'прямоугольный сигнал');
    drawFigure(4, xNum, yGauDFT, yGauDFT_NT, 'ДПФ', 'Гаусс');
end

function y = rectangularPulse(x)
    c = 3;
    y = zeros(size(x));

    y(abs(x) < c) = 1;
end

function y = gaussSignal(x)
    A = 1;
    sigma = 0.5;

    y = A * exp(-(x / sigma).^2);
end

function xArr = dft(x)
    xArr = 0:length(x) - 1;
    temp = -2 * pi * sqrt(-1) * xArr / length(x);
    for i = 1:length(xArr)
		xArr(i) = 0;
		for j = 1:length(x)
			xArr(i) = xArr(i) + x(j) * exp(temp(i) * j);
        end
    end
end

function drawFigure(ind, xNum, y, yImproved, nameAlg, namePulse)
	draw(ind, xNum, y, yImproved, @(x) abs(x), nameAlg, namePulse);
end

function draw(ind, xNum, y, yImproved, f, nameAlg, namePulse)
	subplot(2, 2, ind);
	plot(xNum, f(y) / length(xNum), xNum, f(yImproved) / length(xNum));
	title([nameAlg, ' ', namePulse]);
	legend('С эффектом близнецов', 'Без эффекта близнецов');
end    