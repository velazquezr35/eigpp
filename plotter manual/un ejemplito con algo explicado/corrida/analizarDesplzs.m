

close all
clear variables
clc


% lee los desplazamientos y giros, uz y phiy, de dosarchivos independientes
% los grafica
% calcula la FFT y la grafica


tini = 0;
tfin = 2000;

% ABRIR ARCHIVOS

nEstrats = 0;

correr=false;
if(nEstrats==0)
  str='(echo snl.@1 1 0 && echo d && echo 200035 3 uz.dat && echo 200035 5 phiy.dat && echo s) | curvas';
  [status,cmdout]=dos(str);
  
  [fidUz, msgUz] = fopen('uz.dat','r');
  [fidPhiy, msgPhiy] = fopen('phiy.dat','r');
  correr=isempty(msgUz) && isempty(msgPhiy);
else
  for i=1:nEstrats
    str1=['(echo snl.@' num2str(i) ' 1 0 && echo d &&'];
    str2=[' echo 200035 3 uz' num2str(i) '.dat &&'];
    str3=[' echo 200035 5 phiy' num2str(i) '.dat && echo s) | curvas'];
    [status,cmdout]=dos([str1 str2 str3]);
    
    [aux, msgUz] = fopen(['uz' num2str(i) '.dat'],'r');
    fidUz(i) = aux;
    [aux, msgPhiy] = fopen(['phiy' num2str(i) '.dat'],'r');
    fidPhiy(i) = aux;
    correr(i)=isempty(msgUz) && isempty(msgPhiy);
  end%for
end%if


if( correr )
  
  
  %%%%%%%%%%%%%%%%
  %  LEER DATOS  %
  %%%%%%%%%%%%%%%%
  if(nEstrats==0)
    uzSimpact = fscanf(fidUz, '%e %e', [2 inf]);
    phiySimpact = fscanf(fidPhiy, '%e %e', [2 inf]);
    fclose(fidUz);
    fclose(fidPhiy);
  else
    uzSimpact = [];
    phiySimpact = [];
    for i=1:nEstrats
      uzSimpact = [uzSimpact fscanf(fidUz(i), '%e %e', [2 inf])];
      phiySimpact = [phiySimpact fscanf(fidPhiy(i), '%e %e', [2 inf])];
      fclose(fidUz(i));
      fclose(fidPhiy(i));
    end%for
  end%if
  phiySimpact(2,:) = -rad2deg( phiySimpact(2,:)-phiySimpact(2,1) );
  % NOTAR QUE CAMBIO LA REFERENCIA para considerar solo cambios de ángulo
  % NOTAR   además, que cambio el signo
  %         esto es para graficar correctamente el cambio de ángulo de ataque
  %         al analizar flutter de la pala SNL100-00 (ver notas en reporte de M. Ramis)
  
  
  %%%%%%%%%%%%%%%%%%%%
  %  GRAFICAR DATOS  %
  %%%%%%%%%%%%%%%%%%%%
  fhCheqear = figure();
    hold on
      phUzSimpact = plot(   uzSimpact(1,:),   uzSimpact(2,:), 'Color','k');%,'Marker','o');
      phPhSimpact = plot( phiySimpact(1,:), phiySimpact(2,:), 'Color','r');%,'Marker','o');
    hold off
    grid on
    ahChequear = gca;
    %axis([0 200 -6 12]);
    title('CHEQUEAR RESULTADOS');
    legend( [phUzSimpact  phPhSimpact], ...
              'delpz. vertical', 'giro [°]');
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  CALCULAR Y GRAFICAR FFT  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % uz
  fprintf(['DESPLAZAMIENTO VERTICAL\n']);
  % elimino el transitorio
  ind = find(uzSimpact(1,:)<tini);
  uzSimpact(:,ind)=[];
  ind = find(uzSimpact(1,:)>tfin);
  uzSimpact(:,ind)=[];
  %
  fprintf(['contando picos en la respuesta\n']);
  [pksu,locsu] = findpeaks(uzSimpact(2,:));
  if(~isempty(pksu))
    frequ = (length(pksu)-1) / ( uzSimpact(1,locsu(end)) - uzSimpact(1,locsu(1)) );
    fprintf(['    om_u = ' num2str(frequ) ' Hz = ' num2str(frequ*2*pi) ' rad/s\n']);
  else
    fprintf(['no se encontraron máximos locales en los desplazamientos\n']);
  end%if
  %
  fprintf(['analizando la transformada de Fourier\n']);
  dfu=1/(uzSimpact(1,end)-uzSimpact(1,1));
  frequ=0:dfu:dfu*(length(uzSimpact(1,:))-1);
  tfu=fft(uzSimpact(2,:));
  tamu=ceil(length(tfu)/2);
  tfu=tfu/max(abs(tfu(2:tamu)));
  
  fhu = figure();
    phu = plot(frequ(2:tamu)*2*pi,abs(tfu(2:tamu)),'-ok');
    grid on
    title('FFT para el deslplazamiento vertical');
  
  [pksu,locsu] = findpeaks(abs(tfu),'npeaks',1,'threshold',0.5);
  if(~isempty(pksu))
    fprintf(['    om_u = ' num2str(frequ(locsu)) ' Hz = ' num2str(frequ(locsu)*2*pi) ' rad/s\n']);
    fprintf(['definición en frecuencias df = ' num2str(dfu) ' Hz = ' num2str(dfu*2*pi) ' rad/s\n\n']);
  else
    fprintf(['no se encontraron máximos locales en la transformada de Fourier\n\n']);
  end%if
  
  % phiy
  fprintf(['CAMBIO DE ÁNGULO DE ATAQUE\n']);
  % elimino el transitorio
  ind = find(phiySimpact(1,:)<tini);
  phiySimpact(:,ind)=[];
  ind = find(phiySimpact(1,:)>tfin);
  phiySimpact(:,ind)=[];
  %
  fprintf(['contando picos en la respuesta\n']);
  [pksp,locsp] = findpeaks(phiySimpact(2,:));
  if(~isempty(pksp))
    freqp = (length(pksp)-1) / ( phiySimpact(1,locsp(end)) - phiySimpact(1,locsp(1)) );
    fprintf(['    om_phi = ' num2str(freqp) ' Hz = ' num2str(freqp*2*pi) ' rad/s\n']);
  else
    fprintf(['no se encontraron máximos locales en los desplazamientos\n']);
  end%if
  %
  fprintf(['analizando la transformada de Fourier\n']);
  dfp=1/(phiySimpact(1,end)-phiySimpact(1,1));
  freqp=0:dfp:dfp*(length(phiySimpact(1,:))-1);
  tfp=fft(phiySimpact(2,:));
  tamp=ceil(length(tfp)/2);
  tfp=tfp/max(abs(tfp(2:tamp)));
  
  fhp = figure();
    php = plot(freqp(2:tamp)*2*pi,abs(tfp(2:tamp)),'-or');
    grid on
    title('FFT para el giro');
  
  [pksp,locsp] = findpeaks(abs(tfp),'npeaks',1,'threshold',0.5);
  if(~isempty(pksp))
    fprintf(['    om_phi = ' num2str(freqp(locsp)) ' Hz = ' num2str(freqp(locsp)*2*pi) ' rad/s\n']);
    fprintf(['definición en frecuencias df = ' num2str(dfp) ' Hz = ' num2str(dfp*2*pi) ' rad/s\n\n']);
  else
    fprintf(['no se encontraron máximos locales en la transformada de Fourier\n\n']);
  end%if
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  DESPLAZAMIENTO BASE Y AMPLITUD DE VIBRACIÓN  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if 0
    % uz
    % first time
    [lMin,~,lMax,~]=localMaxMin( uzSimpact(2,:) );
    % filter high freqs. local min and max
    [~,~,lMax,~]=localMaxMin( lMax ); % max
    [lMin,~,~,~]=localMaxMin( lMin ); % min
    uppLimUz = mean(lMax);
    lowLimUz = mean(lMin);
    midPosUz = mean([uppLimUz lowLimUz]);
    axes(ahChequear);
    hold on
      phCheqUppLimUz = plot( [tini uzSimpact(1,end)], [uppLimUz uppLimUz], 'b');
      phCheqLowLimUz = plot( [tini uzSimpact(1,end)], [lowLimUz lowLimUz], 'b');
      phCheqMidPosUz = plot( [tini uzSimpact(1,end)], [midPosUz midPosUz], 'b');
    hold off
    % phiy
    % first time
    [lMin,~,lMax,~]=localMaxMin( phiySimpact(2,:) );
    % filter high freqs. local min and max
    [~,~,lMax,~]=localMaxMin( lMax ); % max
    [lMin,~,~,~]=localMaxMin( lMin ); % min
    uppLimPhiy = mean(lMax);
    lowLimPhiy = mean(lMin);
    midPosPhiy = mean([uppLimPhiy lowLimPhiy]);
    axes(ahChequear);
    hold on
      phCheqUppLimPhiy = plot( [tini phiySimpact(1,end)], [uppLimPhiy uppLimPhiy], 'b');
      phCheqLowLimPhiy = plot( [tini phiySimpact(1,end)], [lowLimPhiy lowLimPhiy], 'b');
      phCheqMidPosPhiy = plot( [tini phiySimpact(1,end)], [midPosPhiy midPosPhiy], 'b');
    hold off
    
    fprintf('DESPLAZAMIENTO BASE Y AMPLITUD DE VIBRACIÓN\n');
    fprintf('variable           uz           phiy [º]\n');
    fprintf('límite inferior   %12.4e %12.4e\n', lowLimUz, lowLimPhiy);
    fprintf('límite superior   %12.4e %12.4e\n', uppLimUz, uppLimPhiy);
    fprintf('amplitud          %12.4e %12.4e\n', abs(uppLimUz-lowLimUz), abs(uppLimPhiy-lowLimPhiy));
    fprintf('posición media    %12.4e %12.4e\n', midPosUz, midPosPhiy);
    fprintf('\n');
  end%if
  
else
  fprintf('\nHUBO PROBLEMAS AL ABRIR LOS ARCHIVOS DE ENTRADA\n\n');
end%if
