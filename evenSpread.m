function initCords = evenSpread(n, sigmaX, sigmaY)
    initCords = zeros(2, 2^n);

    if mod(n,2) == 0
        m = n/2;
        l = n/2;
    end
    if mod(n,2) == 1
        m = (n+1)/2;
        l = (n-1)/2;
    end

    sigmaXinc = sigmaX / 2 ^ l;                 % Value to increment by (distance between points)
    sigmaYinc = sigmaY / 2 ^ m;
    
    xVal = sigmaXinc;
    yVal = sigmaYinc;
    i = 1;
    while xVal < 2 * sigmaX
        while yVal < 2 * sigmaY
            initCords(1,i) = xVal;
            initCords(2,i) = yVal;
            yVal = yVal + 2 * sigmaYinc;
            i = i + 1;
        end
        yVal = sigmaYinc;
        xVal = xVal + 2 * sigmaXinc;
    end
    initCords = initCords - [sigmaX; sigmaY];   % This will center it at 0
end