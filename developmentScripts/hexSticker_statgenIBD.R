# note: run usethis::use_logo() to generate the generate correct size

library(hexSticker)
library(ggplot2)
library(LMMsolver)
library(spam)
library(fields)
library(maps)
library(dplyr)
library(gridExtra)

dat1 <- data.frame(x=1,y=5,group=1)
dat2 <- data.frame(x=1,y=7, group=1)
dat3 <- data.frame(x=1.10,y=5, group=2)
dat4 <- data.frame(x=1.10,y=7, group=2)

ind1 <- rbind(dat1, dat2, dat3, dat4)

dat5 <- data.frame(x=2,y=5, group=3)
dat6 <- data.frame(x=2,y=7, group=3)
dat7 <- data.frame(x=2.10,y=5, group=4)
dat8 <- data.frame(x=2.10,y=7, group=4)

ind2 <- rbind(dat5, dat6, dat7, dat8)

dat9 <- data.frame(x=1.45,y=0.0, group=5)
dat10 <- data.frame(x=1.45,y=0.4, group=5)
dat11 <- data.frame(x=1.55,y=0.0, group=6)
dat12 <- data.frame(x=1.55,y=0.9, group=6)
dat13 <- data.frame(x=1.45,y=0.4, group=7)
dat14 <- data.frame(x=1.45,y=2.0, group=7)
dat15 <- data.frame(x=1.55,y=0.9, group=8)
dat16 <- data.frame(x=1.55,y=2.0, group=8)

ind3 <- rbind(dat9, dat10, dat11, dat12, dat13, dat14, dat15, dat16)
dat <- rbind(ind1, ind2, ind3)

dat$group = as.factor(dat$group)
centromeres <- data.frame(x=c(1.05,2.05,1.5),y=c(6.5,6.5, 1.5))
df_dashed_line <- data.frame(x=c(1.50,1.50),y=c(2.50,5.00))
p <- ggplot(dat,aes(x=x,y=y)) + geom_line(aes(colour=group),lwd=2.0) +
  geom_line(data=df_dashed_line, aes(x=x,y=y), lwd=0.5, lty=3) +
  scale_colour_manual(values=c("green","red","orange","blue","green","orange","red",
                               "blue")) +
  annotate(
    "text", label = "x",
    x = 1.5, y = 6, size = 12, colour = "black"
  ) +
  geom_point(data=centromeres, aes(x=x,y=y),size=1) + xlim(0.9,2.2) + ylim(0.0,7.0) +
 theme_void() + theme_transparent() + theme(legend.position="none")
p

sticker(p, package="statgenIBD", p_size=20, s_x=1.05, s_y=.70, s_width=1.2, s_height=0.9,
        h_fill = "lightyellow", p_color="blue",h_color="green",
        filename="statgenIBD_hexSticker.png")

