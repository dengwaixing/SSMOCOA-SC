function IDX = spectral(W, num_clusters)
% 构建拉普拉斯矩阵
D = diag(sum(W));
L = D - W;
% 特征值分解
[V, ~] = eig(L);
% 选择前k个特征向量

V = V(:,1:num_clusters);
% 对特征向量进行归一化处理
V = normr(V);
% 对低维表达进行聚类
IDX = kmeans(V,num_clusters);
end
