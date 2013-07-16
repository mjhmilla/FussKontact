function contactInfo = calcContactForcePosition(t,x,vParams,vInputs)

tmp = xDotMex(t, x, vParams, vInputs);
contactInfo = tmp((length(x)+1):1:length(tmp));