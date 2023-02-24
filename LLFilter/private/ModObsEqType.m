function [Prediction, Filter, ModEq, ObsEq] = ModObsEqType(FileNames,n_v_e)

% Defines type of Model and Observation Equations

ModEq=FileNames{1};
ObsEq=FileNames{2};

OdeKrylovBarrier=40;
SdeKrylovBarrier=6;

switch nargout(ModEq)
    case 2
       if n_v_e < OdeKrylovBarrier
         Prediction=@Prediction_ODE;
       else
         clear Prediction_ODE_K
         Prediction=@Prediction_ODE_K;
       end
    case 3
       if n_v_e < SdeKrylovBarrier
         Prediction=@Prediction_Auto_Add;
       else
         clear Prediction_Auto_Add_K
         Prediction=@Prediction_Auto_Add_K;
       end         
    case 4
       if n_v_e < SdeKrylovBarrier      
         Prediction=@Prediction_Auto_Mult;
       else
         clear Prediction_Auto_Mult_K
         Prediction=@Prediction_Auto_Mult_K;
       end         
    case 5
       if n_v_e < SdeKrylovBarrier       
         Prediction=@Prediction_NoAuto_Add; 
       else
         clear Prediction_NoAuto_Add_K
         Prediction=@Prediction_NoAuto_Add_K;
       end         
    case 6
       if n_v_e < SdeKrylovBarrier
         Prediction=@Prediction_NoAuto_Mult;
       else
         clear Prediction_NoAuto_Mult_K
         Prediction=@Prediction_NoAuto_Mult_K;
       end         
    otherwise
       disp('Error in Model definition')
end

switch nargout(ObsEq)
    case 3
        Filter=@Filter_Add;
    case 6
        Filter=@Filter_Mult;
    otherwise
        disp('Error in Observation definition')
end

