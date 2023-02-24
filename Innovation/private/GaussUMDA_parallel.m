function [Argument,BestMax] = GaussUMDA_parallel(FitnessFunct,LowerBound,UpperBound,EstimOptions,OptimOptions) 
%   GaussUMDA is an implementation of an estimation of distribution algorithm
%   for the computation of innovation estimators of diffusion processes.
%   This is the parallel implementation of the function GaussUMDA.m
%
% INPUT
% FitnessFunct: Fitness Function, its argument must be a vector of NumbVar variables.
% LowerBound: Vector with the lower bound for the parameters
% UpperBound: Vector with the upper bound for the parameters
% OptimOptions.NumGen: Maximum number of generations in the algorithm
% OptimOptions.PopSizeRate:  PopSizeRate * length(LowerBound) = Population Size
% OptimOptions.TruncValue: Truncation value (cuando T=0, entonces es seleccion proporcional)
%                          para TruncationValue=1 se escoge la poblacion  completa
% OptimOptions.MaxFitness: Maximum permissible value of the fitness function, useful as stop condition  
%
% OUTPUT
% BestMax: Maximum value of the fitness function
% Argument: Argument of the fitness function at the maximum value
%

%  Reference
%   - "Estimation of distribution algorithms for the computation 
%      of innovation estimators of diffusion processes",  
%      González Z., Jimenez J.C., Lozada-Chang L.V., Santana R.
%      Mathematics and Computers in Simulation, 187 (2021) 449–467.

% Copyright (c) González-Arenas Z., 2020.
% Maintainer : zochilita@gmail.com
% Global Optimization UMDA 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, and to permit persons to 
% whom the Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

% k: Generacion en que fue encontrado el maximo de conocerse con anterioridad su valor
% BestArguments: Matriz con la mejor solucion en cada generacion
% FunGen: valor de la funcion objetivo evaluado en la media

rand('seed',sum(clock));

NumVar = length(LowerBound);
PopSize = OptimOptions.PopSizeRate * NumVar;
CantGen = OptimOptions.NumGen;
TruncationValue = OptimOptions.TruncValue;
MaxFitness = OptimOptions.MaxFitness;

k=1;
if TruncationValue < .1         %garantizando un umbral de truncamiento por encima de 0.1
    TruncationValue = .1;
end
if PopSize < 10                 %garantizando un tamanho de poblacion minimo util
    PopSize = 10;
end

%Generando poblacion inicial en la que no haya individuos que den NaN al
%evaluar en la funcion objetivo

FunVal=zeros(PopSize,1);
NewPop=zeros(PopSize,NumVar);
parfor i=1:PopSize
    FVal=NaN;
    while isnan(FVal) || ~isreal(FVal)
      indiv = LowerBound + (rand(1,NumVar).*(UpperBound-LowerBound));
      FVal  = -FitnessFunct(indiv,EstimOptions);
    end
    NewPop(i,:) = indiv;
    FunVal(i) = FVal;
end

disp('     ')
while( (k==1) | (k<=CantGen & Max(k-1)<MaxFitness) ) 
    
       %disp(['k = ' num2str(k)])
       Pop = NewPop;
             
       [Val,Ind] = sort(FunVal);                   % Ordenamiento para el truncation (de menor a mayor)
  
       NanNumb = sum(isnan(Val));                  % cantidad de valores NaN en FunVal
       p = PopSize - NanNumb;                      % cantidad de individuos utiles en la nueva poblacion
       trunc = floor(PopSize*TruncationValue);     % tamanho de la poblacion para aprender el modelo
       
     % Seleccionando la muestra para aprender el modelo probabilistico  
       if p < trunc             
           SelNewPop = Pop(Ind(1:p),:);
           NewFunVal = Val(1:p);
           parfor i=p+1:trunc-1
              NewFVal=NaN;
              NewIndiv=UpperBound;
              while isnan(NewFVal) || ~isreal(NewFVal)
                for j = 1:NumVar
                   mtrunc = vars_mean(j);
                   strunc = vars_sigmas(j);
                   NewIndiv(j) = norminv(rand*(normcdf(UpperBound(j),mtrunc,strunc)-normcdf(LowerBound(j),mtrunc,strunc)) ...
                               + normcdf(LowerBound(j),mtrunc,strunc),mtrunc,strunc);
                end
                NewFVal = -FitnessFunct(NewIndiv,EstimOptions);
              end
              SelNewPop(i,:) = NewIndiv;
              NewFunVal(i) = NewFVal;
           end
           [Val,Ind] = sort(NewFunVal);
           SelPop = SelNewPop(Ind,:);
       else
           SelPop = Pop(Ind(p-trunc+1:p),:); 
       end
                 
       Max(k) = Val(p);                          % Maximo valor en la poblacion
       BestArguments(k,:) = Pop(Ind(p),:);       % Mejor solucion en la poblacion
%        BestArguments(k,:) = Pop(Ind(PopSize),:);       % Mejor solucion en la poblacion
       
     % la poblacion util es de tamanho p         
     % guardando el promedio de ajuste de la poblacion en cada generacion,
       
       cleanPop = Val(1:p);                  % vector de evaluacion en la poblacion descartando los valores NaN
       MeanFitPop(k) = mean(cleanPop);       % promedio de ajuste para la generacion k
       
             
     % Se Calcula el modelo probabilistico
       vars_mean   = mean(SelPop,1);              % media
       vars_cov    = cov(SelPop);                 % matriz de covarianzas 
       vars_sigmas = sqrt(diag(vars_cov))';
       MediaGen(k,:) = vars_mean;               %media de SelPop usada para calcular el modelo probabilistico 
           
       if ( k ~= CantGen)
     % Se genera la proxima poblacion muestreando a partir del modelo 
     % y de acuerdo a su complejidad
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Imponiendo elitismo en la poblacion, se mantiene el 5% mejor de la
     % poblacion anterior en la nueva.
           Elit = floor(p*.05);
           NewPop(1:Elit,:) = Pop(Ind(p - Elit + 1:p),:);
           FunVal(1:Elit) = Val(p - Elit + 1:p);
     % se incluye la media de cada parametro como un individuo en la nueva
     % poblacion.
           NewPop(Elit + 1,:) = vars_mean;
           FVal = -FitnessFunct(vars_mean,EstimOptions);
           if isnan(FVal) || ~isreal(FVal) 
              FVal=-inf;
           end
           FunVal(Elit + 1) = FVal; 
           
       % debemos sacar la media por generacion (ya esta guardada) y su funval, porque martica
       % lo quiere ver, comparar FunGen(k) con Max(k)
       %    FunGen(k)    = FunVal(Elit + 1);
           
     % Generando el resto de la nueva poblacion
           parfor i = Elit + 2: PopSize
           % Imponiendo que el indiv generado este dentro de los
           % intervalos segun una normal truncada
              FVal = NaN;
              indiv=UpperBound;
              while isnan(FVal) || ~isreal(FVal) 
                 for j=1:NumVar
                    mtrunc = vars_mean(j);
                    strunc = vars_sigmas(j);
                    indiv(j) = norminv(rand*(normcdf(UpperBound(j),mtrunc,strunc) - normcdf(LowerBound(j),mtrunc,strunc)) ...
                             + normcdf(LowerBound(j),mtrunc,strunc),mtrunc,strunc);
                 end
                 FVal = -FitnessFunct(indiv,EstimOptions);
               end
              NewPop(i,:) = indiv;
              FunVal(i) = FVal;
           end
       end
       
       disp([' Generation ',num2str(k),'   Fitness Function = ',num2str(-Max(k))])
      
       k=k+1; 
 end

[BestMax,I] = max(Max);
Argument = BestArguments(I,:);

return 


    
