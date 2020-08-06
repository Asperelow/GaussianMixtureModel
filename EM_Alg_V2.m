clear all
close all

randNumsPerCluster = 10000;                     % Number of random points per cluster
n = 4;
k = 2^n;                                         % Number of centroids

           
randPoints = randn([2 randNumsPerCluster]);     % Creates normally distributed clusters

mag1 = [2 0;0 1];
theta1 = 45 * (pi/180); 
angle1 = [cos(theta1) sin(theta1);-sin(theta1) cos(theta1)];
cluster1 = (angle1 * mag1 * randPoints) + [5; 5];
mag2 = [1 0;0 2];
theta2 = -30 * (pi/180);
angle2 = [cos(theta2) sin(theta2);-sin(theta2) cos(theta2)];
cluster2 = (angle2 * mag2 * randPoints) + [10;10];
mag3 = [3 0;0 0.8];
theta3 = -70 * (pi/180);
angle3 = [cos(theta3) sin(theta3);-sin(theta3) cos(theta3)];
cluster3 = (angle3 * mag3 * randPoints) + [0;-5];
mag4 = [3 0;0 0.8];
theta4 = -70 * (pi/180);
angle4 = [cos(theta3) sin(theta3);-sin(theta3) cos(theta3)];
cluster4 = (angle3 * mag3 * randPoints) + [3;-6];
cluster5 = (randPoints) + [8;-3.5];
randNumsMat =[cluster1, cluster2, cluster3 ,cluster4,cluster5];
randNums = length(randNumsMat(2));   

figure(1)
scatter(randNumsMat(1,:), randNumsMat(2,:),'m.')

randPointsx = randNumsMat(1,:);
randPointsy = randNumsMat(2,:);
clear mag1 mag2 mag3 angle1 angle2 angle3 randNumsPerCluster
clear theta1 theta2 theta3

covData = cov(randNumsMat');
sigmaX = covData(1,1);
sigmaY = covData(2,2);
mu = evenSpread(n, sigmaX, sigmaY);
mu = mu + [mean(randPointsx); mean(randPointsy)]    % evenSpread returns values centered at zero, they need to be shifted to
mu_x = mu(1,:);                                     % the mean of the input data
mu_y = mu(2,:);                 
sigma{k} = 0;
for i = 1:k
    sigma{i} = covData;
end

goodFit = false;
iteration = 1;
fit = 0;
while goodFit == false
                     
    for i = 1:k
        probMat(i,:) = mvnpdf(randNumsMat',[mu_x(i) mu_y(i)],sigma{i});     % Probability for any point existing on any gaussian
    end
    weightMat = probMat(:,:) ./ sum(probMat(:,:));                          % Probability of any point belonging to any centroid compared to others
    mu = [mu_x(:) mu_y(:)];
    oldMu = mu;
    for i = 1:k
        mu_x(i) = sum(weightMat(i,:) .* randNumsMat(1,:))/ sum(weightMat(i,:));   % Updates x and y vectors for mu
        mu_y(i) = sum(weightMat(i,:) .* randNumsMat(2,:))/ sum(weightMat(i,:));
    
        xy_cov(i) = sum(weightMat(i,:) .* ((randPointsx - mu_x(i)) .* (randPointsy - mu_y(i))));
        x_var(i) = sum(weightMat(i,:) .* (randPointsx - mu_x(i)) .^ 2);
        y_var(i) = sum(weightMat(i,:) .* (randPointsy - mu_y(i)) .^ 2);
        sigma{i} = [x_var(i) xy_cov(i); xy_cov(i) y_var(i)] ./ sum(weightMat(i,:));
    end
    
    gaussHeight = sum(weightMat,2) / sum(sum(weightMat));       % Assigns a height to each gaussian
    likelyhood = sum(probMat');
    fit =  [fit, sum(gaussHeight(:) .* likelyhood(:))];
    
    if iteration > 1
        if fit(iteration) - fit(iteration-1) < 0.005
            goodFit = true;
        end
    end
    x1 = -20:0.3:20;
    x2 = -20:0.3:20;
    [X1,X2] = meshgrid(x1,x2);
    X = [X1(:) X2(:)];
    y = 0;
    for i = 1:k
        y = y + mvnpdf(X,[mu_x(i) mu_y(i)],sigma{i}) * gaussHeight(i);
    end
    y = reshape(y,length(x2),length(x1));
    figure(2)
    surf(x1,x2,y)
    pause(0.1)
    mu = [mu_x(:) mu_y(:)];
    iteration = iteration + 1;
end
figure(3)
t = 1:iteration;
plot(t,fit(:))
