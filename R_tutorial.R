
######### 1 Vectors and assignment
x <- c(10.4, 5.6, 3.1, 6.4, 21.7)
assign("x", c(10.4, 5.6, 3.1, 6.4, 21.7))
c(10.4, 5.6, 3.1, 6.4, 21.7) -> x
y <- c(x, 0, x)

#Vector arithmetic
#Mathematical operators : +, -, *, / and ^
v <- 2*x + y + 1

#Mathematical function: mean(), length(), sum(), min(), max(), log()
sum((x-mean(x))^2)/(length(x)-1)

#Generating vector sequencing
1:10
seq(length=50, from=0, to=20)
seq(length=50, from=-5, by=0.2)
rep(x, times=5)
rep(x, each=5)

#Logical vectors
#Logical operators <, <=, >, >=, ==, !=, & (and), | (or)

vlogical <- x > 13
(x<5) | (x>20)

#Missing Values
z <- c(1:3,NA)
is.na(z)


#Character Vectors
cv <- c("X1","Y1","X2","Y2","X3","Y3")
cv2 <- paste("X", 1:10, sep="")

#Subsetting a vector
cv2[3]
cv2[c(1,10)]
cv2[-4]
cv2[cv2 == "X1"]

######## 2 Matrix
set.seed(2342)
v1 <- sample(1:10, 10, replace = TRUE)
m1 <- matrix(v1, ncol=2, nrow=5)

#Matrix Operations
m1 * 2  #Multiplication by a constant
m1 * m1 
m1 ^ 2

v2 <- sample(1:10, 5, replace = TRUE)
m1 * v2

#Adding a column 
m2 <- cbind(m1,v2)
colnames(m2) <- c("Col1", "Col2", "Col3")
dim(m2)

#Adding a row 
v3 <- c(4,6,1)
m3 <- rbind(m2, v3)
rownames(m3) <- 1:nrow(m3)
rownames(m3) <- NULL

#Subset a matrix
m3[1,]
m3[,2]
m3[2,2]

######## 3 Data frame
v1 <- seq(0,10,2)
v2 <- rep(c("A","B"), 3)
v3 <- rep(c(TRUE,FALSE), each=3)

df <- data.frame(v1,v2,v3)
colnames(df) <- c("Numeri","Caratteri","Logici")

#Subset a data frame 
df$Numeri
df[,1]
df[df$Logici == TRUE, ]
df[df$Caratteri == 'A', ]

#Add a column 
v4 <- c("bianco","grigio chiaro","grigio","grigio scuro","grigio molto scuro","nero")
df$Colori <- v4

#Ordinal variable become factors
df$Colori <- factor(df$Colori, levels=c("bianco","grigio chiaro","grigio","grigio scuro", "grigio molto scuro","nero"))






