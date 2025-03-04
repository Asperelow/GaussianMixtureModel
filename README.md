# GaussianMixtureModel

### Expectation Maximization (EM) Algorithm

The EM algorithm is  more advanced way to classify groups of datapoints then an algorithm like k-means. K-means will return the average value of all points in an expected cluster, whereas the EM algorithm assigns Gaussian distributions to each cluster, allowing for a more accurate classification of groups (in case groups fit in more oval-like shapes, not circles.) The EM algorithm can be inialized in a number of ways, like using k-means to chose the initial groups, randomly placing centroid over datapoints like you would k-means (though, this has its own problems), and casting a "net" of points eveny spaced with a covariance matrix equal to the covariance of every datapoint. After the algorithm is initiallized, every point is compared to every centroid to determine which centroid is most likely to belong to any point. The likelyhood of any given point belonging to any given centroid is used to update the centroids. This is repeated until the centroids converge.

#### EM_Alg_V2

In the attached code, first the number of randomally generated datapoints are declared, along with the number of centroids (2^n). After that, the random points are created and the "net" of initial sampled points are returned from the evenSpread function. The input to the evenSpread function include the number of points in the height or width, the variance of X, and the Variance of Y. Then,the loop will iterate through the above mentioned process for the EM algorithm until it converges, which is determined when the most any centroid moves in the x or y directions is less then 0.005.. Then, it will display how well the GMM fits the data overtime.

Version 2 (V2) is only capable of 2-Dimensional GMMs, V3 is in progress for any number of dimension GMMs. 
