#--------------------------------------------------
# method 1: P = exp(Q t)
#-------------------------------------------------- 
library(rmutil)  # contains mexp()
path.probs1 <- function(t, Q)
{
    P <- mexp(Q * t)
    
    ans <- c(P[1,3], P[1,4])
    names(ans) <- c("Wa", "Wb")

    return(ans)
}

# Note: The matrix exponential function is also available in the msm library as
# MatrixExp(), but that seems to perform poorly for large t.  So use mexp here, 
# but also check with method 2 when assuming irreversibility.

#--------------------------------------------------
# method 2: P_14a and P_14b
#-------------------------------------------------- 
path.probs2 <- function(t, q13, q14, q34)
{
    q134 <- q13 + q14
    e1t <- exp(-q134 * t)

    ppa <- q14 / q134 * (1 - e1t)

    if (abs(q34 - q13 - q14) > 1e-10)
    {
        ppb <- q13 / q134 * 
              (1 - 1/(q134-q34) * (q134 * exp(-q34*t) - q34 * e1t) )
    } else
    {
        ppb <- q13 / q34 * (1 - (1 + q34 * t) * exp(-q34 * t))
    }

    ans <- c(ppa, ppb)
    names(ans) <- c("Wa", "Wb")

    return(ans)
}

#--------------------------------------------------
# usage
#-------------------------------------------------- 

# 1 = SI, D
# 3 = SC, D
# 4 = SC, P
# Wa = 1 -> 4
# Wb = 1 -> 3 -> 4

# MAP values from EstimatesIDCDPnodelta.csv
q14 <- 0.0681          # ID -> CP, rho_I
q13 <- 0.2066          # ID -> CD, q_IC
q34 <- 0.0370          # CD -> CP, rho_C
q41 <- q31 <- q43 <- 0 # irreversible

# Q is only needed for method 1
Q <- matrix(c(
        -q13 - q14,  q13,        q14,  0,
        q31,        -q31 - q34,  0,    q34,
        q41,         0,         -q41,  0,
        0,           q43,        0,   -q43
        ), byrow=TRUE, nrow=4)

tmax <- 50
times <- seq(0, tmax, tmax/1000)

# pp$Wa is W_a(t), and pp$Wb is W_b(t)
# choose one of these:
pp <- as.data.frame(t(sapply(times, path.probs2, q13, q14, q34)))
pp <- as.data.frame(t(sapply(times, path.probs1, Q)))

# plot
plot(times, pp$Wa, type="l", lwd=3, xlab="elapsed time (My)", ylab="probability", ylim=c(0,1))
lines(times, pp$Wb, lwd=3, lty=2)
legend("topleft", c("ID to CP", "ID to CD to CP"), lty=c(1, 2), lwd=3)

plot(times, pp$Wa/pp$Wb, log="y", type="l", lwd=3, xlab="elapsed time (My)", ylab="one-step : two-step", ylim=c(0.3, 100))
abline(h=1, lty=3)

plot(times, pp$Wa / (pp$Wa + pp$Wb), type="l", lwd=3, xlab="elapsed time (My)", ylab="proportion one-step", ylim=c(0,1))
abline(h=0.5, lty=3)

### Try with diversification ###

r1 <-  0.4607
r3 <-  0.0684
r4 <- -0.0877

Qd <- matrix(c(
       -q13 - q14 + r1,  q13,             q14,       0,
        q31,            -q31 - q34 + r3,  0,         q34,
        q41,             0,              -q41 + r4,  0,
        0,               q43,             0,        -q43 + r4
        ), byrow=TRUE, nrow=4)

ppd <- as.data.frame(t(sapply(times, path.probs1, Qd)))

plot(times, ppd$Wa, type="l", lwd=3, xlab="elapsed time", ylab="relative 'probability'", log="y")
lines(times, ppd$Wb, lwd=3, lty=2)
legend("topleft", c("SI/D to SC/P", "SI/D to SC/D to SC/P"), lty=c(1, 2), lwd=3)

plot(times, ppd$Wa/ppd$Wb, log="y", type="l", lwd=3, xlab="elapsed time (My)", ylab="one-step : two-step", ylim=c(0.3, 100))
abline(h=1, lty=3)

plot(times, ppd$Wa / (ppd$Wa + ppd$Wb), type="l", lwd=3, xlab="elapsed time (My)", ylab="proportion one-step", ylim=c(0,1))
abline(h=0.5, lty=3)

#--------------------------------------------------
# Final plot
#--------------------------------------------------

pdf(file="../manuscript/pathways.pdf", width=8, height=7, family="Times")
par(mfrow=c(2,2), mar=c(4,4,0,0.5), oma=c(0,0,2,0), cex=1.05, cex.lab=1.1)

# trans only
plot(times, pp$Wa, type="l", lty=2, lwd=3, xlab="", ylab="contribution of each pathway", ylim=c(0,1))
lines(times, pp$Wb, lwd=3, lty=3)
# legend("topleft", c("one-step: ID to CP", "two-step: ID to CD to CP"), lty=c(2, 3), lwd=3)
legend("topleft", c("one-step", "two-step"), lty=c(2, 3), lwd=3)
mtext("without diversification", line=1, font=3, cex=1.3)

# with div
plot(times, ppd$Wa, type="l", lty=2, lwd=3, xlab="", ylab="", log="y", yaxt="n")
lines(times, ppd$Wb, lwd=3, lty=3)
ticks <- axTicks(2, log=T)
axis(2, at = ticks, labels = formatC(ticks, format = "g"))
# legend("topleft", c("one-step: ID to CP", "two-step: ID to CD to CP"), lty=c(2, 3), lwd=3)
legend("topleft", c("one-step", "two-step"), lty=c(2, 3), lwd=3)
mtext("with diversification", line=1, font=3, cex=1.3)

# trans only
plot(times, pp$Wa / (pp$Wa + pp$Wb), type="l", lwd=3, xlab="elapsed time (My)", ylab="proportion of one-step", ylim=c(0,1))
abline(h=0.5, lwd=2, lty=3)

# with div
plot(times, ppd$Wa / (ppd$Wa + ppd$Wb), type="l", lwd=3, xlab="elapsed time (My)", ylab="", ylim=c(0,1))
abline(h=0.5, lwd=2, lty=3)

dev.off()
