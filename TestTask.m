% close all;
clear all;
clc;

%% Параметры эксперимента
blockLength = 500;
dataBlockCount = 10;
sampleRate = 200e3;
toneFreq = 50e3;
freqError = 10e3; % Увеличили частотную расстройку до 10 КГц
snrDbRange = -20:2:20;
seedCount = 1000;

%% Генерация передаваемого сигнала
% Последовательность для коррекции частоты: комплексная экспонента с частотой toneFreq
freqSeq = exp(1j*2*pi*toneFreq/sampleRate*(1:blockLength));

% Синхронизирующая последовательность - случайный комплексный сигнал нормированный по мощности
% Используется для корреляции
rng(1);
synchSeq = randn(1, blockLength) + 1j*randn(1, blockLength);
synchSeq = synchSeq / sqrt(mean(abs(synchSeq).^2));

% Последовательность данных - модуляция BPSK: 1 = бит 1, -1 = бит 0
dataSeq = randi([0 1], 1, blockLength*dataBlockCount) * 2 - 1;

% Формируем полный сигнал
txSignal = [freqSeq synchSeq dataSeq];

%% Запуск эксперимента
syncFailed = zeros(seedCount, length(snrDbRange));
% Перебираем все значения SNR
for snrDbIdx = 1:length(snrDbRange)
    snrDb = snrDbRange(snrDbIdx);
    
    % Перебираем все значения seed для статистики
    for seedIdx = 1:seedCount
        rng(seedIdx);
        
        %% Прохождение через канал
        % Применяем частотный сдвиг из-за рассогласования несущих
        rxSignal = txSignal .* exp(1j*2*pi*freqError/sampleRate*(1:length(txSignal)));
        
        % Генерируем AWGN шум из распределения CN(0,1)
        noise = (randn(size(rxSignal)) + 1j * randn(size(rxSignal))) * sqrt(1/2);
        
        % Добавляем шум с заданным SNR
        rxSignal = rxSignal + noise * db2mag(-snrDb);
        
        %% Синхронизация
        % *** Грубая временная и частотная синхронизация по freqSeq ***
        
        % Вычисляем разность фаз между последовательными отсчетами
        phaseDiff = rxSignal(2:end) .* conj(rxSignal(1:end-1));
        
        % Применяем скользящее среднее с окном длиной blockLength (длина freqSeq)
        windowSize = blockLength;
        slidingAvg = abs(filter(ones(1, windowSize), windowSize, phaseDiff));
        
        % Дополняем slidingAvg до длины rxSignal
        slidingAvg = [zeros(1,1), slidingAvg];
        
        % Находим пик в slidingAvg, соответствующий freqSeq
        [~, coarseTimeIdx] = max(slidingAvg);
        
        % Оцениваем частотную расстройку по freqSeq
        % Извлекаем участок freqSeq из принятого сигнала
        if coarseTimeIdx + blockLength -1 <= length(rxSignal)
            freqSeqRx = rxSignal(coarseTimeIdx : coarseTimeIdx + blockLength -1);
        else
            % Если выход за пределы массива, помечаем как неуспешную синхронизацию
            syncFailed(seedIdx, snrDbIdx) = 1;
            continue;
        end
        
        % Вычисляем среднюю разность фаз
        phaseDiffFreqSeq = angle(freqSeqRx(2:end) .* conj(freqSeqRx(1:end-1)));
        estFreqError = mean(phaseDiffFreqSeq) * sampleRate / (2*pi);
        
        % Компенсируем частотную расстройку
        freqOffsetCompensation = exp(-1j*2*pi*estFreqError/sampleRate*(1:length(rxSignal)));
        rxSignalComp = rxSignal .* freqOffsetCompensation;
        
        % *** Точная временная синхронизация по synchSeq ***
        
        % Корреляция с синхронизирующей последовательностью
        correlation = conv(rxSignalComp, conj(synchSeq(end:-1:1)));
        
        % Поиск максимума корреляции
        [~, maxIdx] = max(abs(correlation));
        
        % Определение начала синхронизации
        syncStart = maxIdx - length(synchSeq) + 1;
        
        % Проверка точности синхронизации
        refSyncStart = length(freqSeq) + 1;
        syncFailed(seedIdx, snrDbIdx) = abs(syncStart - refSyncStart) > 3;
    end
end

%% Построение графика результатов
fig = figure; semilogy(snrDbRange, mean(syncFailed,1),'Linewidth',2);
xlabel("SNR, dB"); ylabel("Probability of failed sync");
ylim([1e-3 1])
grid on; title("Initial frequency mismatch is "+string(round(freqError/1e3))+" KHz");
saveas(fig, "res_"+string(round(freqError/1e3))+".png");