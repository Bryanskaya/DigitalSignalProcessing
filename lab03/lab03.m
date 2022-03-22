function lab03
	c = 3.0;
	sigma = 1.0;

    n = input('Введите количество точек: ');
    step = input('Введите шаг: ');

    t_max = step * (n - 1) / 2;
	t = -t_max:step:t_max;

	x1 = [rectangularPulse(t, c) zeros(1, length(t))];
	x2 = [gaussSignal(t, sigma) zeros(1, length(t))];
	x3 = [rectangularPulse(t, c/2) zeros(1, length(t))];
	x4 = [gaussSignal(t, sigma/2) zeros(1, length(t))];

	% Свертка
	y1 = ifft(fft(x1).* fft(x2)) * step; % Прямоугольный + Гаусс
	y2 = ifft(fft(x1).* fft(x3)) * step; % Прямоугольный + прямоугольный
	y3 = ifft(fft(x2).* fft(x4)) * step; % Гаусс + Гаусс

	start = fix((length(y1) - length(t)) / 2);
	y1 = y1(start+1:start+length(t));
	y2 = y2(start+1:start+length(t));
	y3 = y3(start+1:start+length(t));

    figure(1);
	graph_figure(1, t, x1, x2, y1, 'Прямоугольный + Гаусс');
	graph_figure(2, t, x1, x3, y2, 'Прямоугольный + прямоугольный');
	graph_figure(3, t, x2, x4, y3, 'Гаусс + Гаусс');
end

function graph_figure(ind, t, x1, x2, y, tit)
    subplot(2, 2, ind);
	plot(t, x1(1:length(t)), t, x2(1:length(t)), t, y);
	title(['Свёртка ', tit]);
	legend('Сигнал 1', 'Сигнал 2', 'Свёртка');
end

function y = rectangularPulse(x, c)
	y = zeros(size(x));
	y(abs(x) - c < 0) = 1;
end

function y = gaussSignal(x, sigma)
    A = 1;
	y = A * exp(-(x / sigma) .^ 2);
end