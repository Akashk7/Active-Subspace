function [srgt_WAS, PRESSRMS_WAS, Y_hat_WAS, eXV_WAS, srgtOPT_WAS, predvar_WAS, R2_pred] = metamodel_WAS(x,y)

X_WAS=x;
Y_WAS=y;
[srgt_KRG, PRESSRMS_KRG, eXV_KRG, srgtOPT_KRG] = build_KRG_SRGT(X_WAS,Y_WAS);
[srgt_PRS, PRESSRMS_PRS, eXV_PRS, srgtOPT_PRS] = build_PRS_SRGT(X_WAS,Y_WAS);
[srgt_RBF, PRESSRMS_RBF, eXV_RBF, srgtOPT_RBF] = build_RBF_SRGT(X_WAS,Y_WAS);
[srgt_WAS, PRESSRMS_WAS, eXV_WAS, srgtOPT_WAS,] = build_WAS_SRGT(X_WAS,...
                                       srgt_KRG, srgtOPT_KRG, eXV_KRG,...
                                       srgt_PRS, srgtOPT_PRS, eXV_PRS,...
                                       srgt_RBF, srgtOPT_RBF, eXV_RBF);
Y_hat_WAS = srgtsWASEvaluate(X_WAS, srgt_WAS);
predvar_WAS = srgtsWASPredictionVariance(X_WAS, srgt_WAS);

%  R^2 predicted
% sum of square about mean
   SST=sum((Y_hat_WAS-mean(Y_hat_WAS)).^2);
%  sum of square due to error
   SSE_pred=((sum(eXV_WAS.^2)));
   R2_pred=1-(SSE_pred/SST);
   