# install.packages("quantmod"); install.packages("lattice")
# install.packages("timeSeries"); install.packages("rugarch")
library(quantmod); library(lattice); library(timeSeries); library(rugarch)

getSymbols("AMZN", src = "google", from="2004-01-01")  #loads data for amazon shares close value
AMZNReturns = diff(log(Cl(AMZN)))  #differenced log of close values
AMZNReturns[1] = 0  #removes the NA value on the first observation

windowLength = 250
foreLength = length(AMZNReturns) - windowLength
forecasts <- vector(mode="character", length=foreLength)  #will contain forecasted values

for (d in 0:foreLength) {   # don't run every time!
    AMZNReturnsOffset = AMZNReturns[(1+d):(windowLength+d)]  
    
    final.aic <- Inf
    final.order <- c(0,0,0)
    
    for (p in 0:5) for (q in 0:5) { #FIX
        arimaFit = tryCatch( arima(AMZNReturnsOffset, order=c(p, 0, q)), 
                             error=function( err ) FALSE, 
                             warning=function( err ) FALSE )
        
        if( !is.logical( arimaFit ) ) {
            current.aic <- AIC(arimaFit)
            if (current.aic < final.aic) {
                final.aic <- current.aic
                final.order <- c(p, 0, q)
                final.arima <- arima(AMZNReturnsOffset, order=final.order)
            }
        } else {
            next
        }
    }
    
    spec = tryCatch(ugarchspec(
        variance.model=list(garchOrder=c(1,1)),
        mean.model=list(armaOrder=c(final.order[1], final.order[3]), include.mean=T),
        distribution.model="sged"),
        error=function(e) e, warning=function(w) w)
 
    fit = tryCatch(
        ugarchfit(
            spec, AMZNReturnsOffset, solver = 'hybrid'
        ), error=function(e) e, warning=function(w) w
    )
    
    if(is(fit, "warning")) {
        forecasts[d+1] = paste(index(AMZNReturnsOffset[windowLength]), 1, sep=",")
        print(paste(index(AMZNReturnsOffset[windowLength]), 1, sep=","))
    } else {
        fore = ugarchforecast(fit, n.ahead=1)
        ind = fore@forecast$seriesFor
        forecasts[d+1] = paste(colnames(ind), ifelse(ind[1] < 0, -1, 1), sep=",")
        print(paste(colnames(ind), ifelse(ind[1] < 0, -1, 1), sep=",")) 
    }
}

library(reshape2)
strategy_df <- colsplit(string = forecasts, pattern = ",", names=c("date", "strategy"))
strategy_df_2 <- data.frame(date = strategy_df$date[-1], 
                            strategy = strategy_df$strategy[1:nrow(strategy_df)-1])
forecasts_new <- paste(strategy_df_2$date, strategy_df_2$strategy, sep = ",")
write.table(forecasts_new, file="forecasts_AMZN_WL250p5q5.csv", sep = ",",
            row.names = F, col.names = F, quote = F)

AMZNArimaGarch = as.xts(
    read.zoo(file = "forecasts_AMZN_WL250p5q5.csv", format = "%Y-%m-%d", header = F, sep = ",")
)

AMZNIntersect = merge(AMZNArimaGarch, AMZNReturns, all=F )
AMZNArimaGarchReturns = AMZNIntersect[,1] * AMZNIntersect[,2]

AMZNArimaGarchCurve = log( cumprod( 1 + AMZNArimaGarchReturns ) )
AMZNBuyHoldCurve = log( cumprod( 1 + AMZNIntersect[,2] ) )
AMZNCombinedCurve = merge( AMZNArimaGarchCurve, AMZNBuyHoldCurve, all=F )

returns_plot <- xyplot( 
    AMZNCombinedCurve,
    superpose=T,
    col=c("darkred", "darkblue"),
    lwd=2,
    key=list( 
        text=list(
            c("ARIMA+GARCH", "Buy & Hold")
        ),
        lines=list(
            lwd=2, col=c("darkred", "darkblue")
        )
    )
)
trellis.device(device = "png", filename = "Returns_AMZN_WL250p5q5.png")
print(returns_plot)
dev.off()

print(returns_plot)
