
library(readr)
model_file_name <- "LorenzGenerate.bi"
writeLines(read_file(model_file_name))

library('rbi')
library(ggplot2)

Lorenz <- bi_model(model_file_name)

T <- 10.0
nObs <- 100
init_parameters <- list(X = 1, Y = 1, Z = 1)

synthetic_dataset <- bi_generate_dataset(end_time=T, model=Lorenz,
                                         init=init_parameters,
                                         noutputs = nObs)

synthetic_data <- bi_read(synthetic_dataset)
synthetic_df <- as.data.frame(synthetic_data)
tail(synthetic_df)

p <- ggplot(synthetic_df, aes(X.time)) +
    geom_path(aes(y = X.value, colour="alpha 16.0")) +
    theme(legend.position="bottom") +
    ggtitle("Lorenz") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Time") +
    ylab("X Value")
ggsave(filename = "diagrams/xpath.svg", plot = p)

path0 <- ggplot() +
    theme(legend.position="bottom") +
    ggtitle("Lorenz") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Time") +
    ylab("Value")


set.seed(42)

T <- 20.0

for (i in c("red", "blue", "green")) {
    init_parameters <- list(X = 1 + rnorm(1,0.0,0.01),
                            Y = 1 + rnorm(1,0.0,0.01),
                            Z = 1 + rnorm(1,0.0,0.01))

    synthetic_dataset <- bi_generate_dataset(end_time=T, model=Lorenz,
                                             init=init_parameters,
                                             noutputs = nObs)

    synthetic_data <- bi_read(synthetic_dataset)
    synthetic_df <- as.data.frame(synthetic_data)

    path0 <- path0 +
        geom_line(data = synthetic_df, aes(x = X.time, y = X.value), color = i)
}

ggsave(filename = "diagrams/xpath4.svg", plot = path0)

model_file_name <- "LorenzState.bi"
writeLines(read_file(model_file_name))

LorenzState <- bi_model(model_file_name)

bi_state_model <- libbi(model=LorenzState)
bi_state <- filter(bi_state_model,
                   nparticles = 8192,
                   nthreads = 1,
                   end_time = T,
                   obs = synthetic_dataset,
                   init = init_parameters,
                   ess_rel = 1,
                   sample_obs = TRUE)

bi_file_summary(bi_state$output_file_name)
bi_state
summary(bi_state)

output <- bi_read(bi_state)
logw <- xtabs(value ~ time + np, data = output$logweight, addNA = TRUE)
X <- output$X$value
Y <- output$Y$value
Z <- output$Z$value
A <- output$ln_alpha$value

log2normw <- function(lw){
  w <- exp(lw - max(lw))
  return(w / sum(w))
}

w = t(apply(X=logw, MARGIN=1, FUN=log2normw))
Xmeans = apply(X = X*w, MARGIN=1, FUN=sum)
Ymeans = apply(X = X*w, MARGIN=1, FUN=sum)
Zmeans = apply(X = Z*w, MARGIN=1, FUN=sum)
Ameans = apply(X = A*w, MARGIN=1, FUN=sum)


synthetic_data <- bi_read(synthetic_dataset)
X_original <- synthetic_data$X$value
Y_original <- synthetic_data$Y$value
Z_original <- synthetic_data$Z$value


synthetic_df <- as.data.frame(synthetic_data)
synthetic_df$Xmeans <- Xmeans
synthetic_df$Ymeans <- Ymeans
synthetic_df$Zmeans <- Zmeans
synthetic_df$Ameans <- Ameans

pAmeans <- ggplot(synthetic_df, aes(X.time)) +
           geom_path(aes(y = exp(Ameans), colour="Ameans")) +
           theme(legend.position="bottom") +
           ggtitle("Lorenz") +
           theme(plot.title = element_text(hjust = 0.5)) +
           ylim(0.0, max(exp(Ameans))) +
           xlab("Time") +
           ylab("Value")


ggsave(filename = "diagrams/xpath5.svg", plot = pAmeans)

dataset_list <- list()
parameters_list <- list()

for (i in c(1,2,3)) {
    init_parameters <- list(X = 1 + rnorm(1,0.0,0.01),
                            Y = 1 + rnorm(1,0.0,0.01),
                            Z = 1 + rnorm(1,0.0,0.01))

    parameters_list[[i]] <- init_parameters
    synthetic_dataset <- bi_generate_dataset(end_time=T, model=Lorenz,
                                             init=init_parameters,
                                             noutputs = nObs)

    dataset_list[[i]] <- synthetic_dataset
}

X_list <- list()
Y_list <- list()
Z_list <- list()
A_list <- list()

 for (i in c(1,2,3)) {
    bi_state <- filter(bi_state_model, nparticles = 8192, nthreads = 1, end_time = T, obs = dataset_list[[i]], init = parameters_list[[i]], ess_rel = 1, sample_obs = TRUE)
    output <- bi_read(bi_state)
    logw <- xtabs(value ~ time + np, data = output$logweight, addNA = TRUE)
    w = t(apply(X=logw, MARGIN=1, FUN=log2normw))
    X <- output$X$value
    Y <- output$Y$value
    Z <- output$Z$value
    A <- output$ln_alpha$value
    X_list[[i]] = apply(X = X*w, MARGIN=1, FUN=sum)
    Y_list[[i]] = apply(X = X*w, MARGIN=1, FUN=sum)
    Z_list[[i]] = apply(X = Z*w, MARGIN=1, FUN=sum)
    A_list[[i]] = apply(X = A*w, MARGIN=1, FUN=sum)
}


path2 <- ggplot() +
    theme(legend.position="bottom") +
    ggtitle("Lorenz") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Time") +
    ylab("Value")


for (i in c(1,2,3)) {
        synthetic_data <- bi_read(dataset_list[[i]])
        synthetic_df <- as.data.frame(synthetic_data)
        synthetic_df$Ameans <- exp(A_list[[i]])
        path2 <- path2 + geom_line(data = synthetic_df,
                                   aes(x = X.time, y = Ameans), color = "blue")
}

ggsave(filename = "diagrams/xpath7.svg", plot = path2)

x <- list()

for (i in c(1:3)) {
        x[[i]] <- tail(exp(A_list[[i]]), n = 50)
}

for (i in 1:3) print(mean(x[[i]]))

for (i in 1:3) print(sd(x[[i]]))
