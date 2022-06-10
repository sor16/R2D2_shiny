#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(extraDistr)
library(ggplot2)
library(shinydashboard)
library(dplyr)

gg_theme <- function(){
    theme_classic()+
    theme(
        plot.title=element_text(size=18,face="bold"),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16)
    )
}

gg_histogram_style <- function(p,title,n_bins){
    p + geom_histogram(color="#000000",fill="#0099F8",bins=n_bins) +
        scale_y_continuous(expand=c(0, 0)) +
        labs(title=title,x="x",y="Count") +
        gg_theme()
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    withMathJax(),
    # section below allows in-line LaTeX via $ in mathjax.
    # tags$div(HTML("
    #             MathJax.Hub.Config({
    #             tex2jax: {inlineMath: [['$','$']}
    #             });
    #             ")),
    # Application title
    titlePanel("R2-D2 prior distribution"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("p",
                        "Number of regression parameters",
                        value=4, 
                        min = 1,max=100,step=1),
            sliderInput("a",
                        "$R^2$ shape parameter - a",
                        value =1, 
                        min = 0,max=10,step=0.1),
            sliderInput("b",
                         "R^2 shape parameter - b",
                         value =1, 
                         min = 0,max=10,step=0.1),
            sliderInput("a_pi",
                         "Concentration parameter - a_pi",
                         value = 1, 
                         min = 0,max=10,step=0.1),
            sliderInput("N",
                        "Number of samples",
                        value = 1000, 
                        min = 0,max=10000),
            sliderInput("bins",
                        "Number of bins:",
                        min = 0,
                        max = 100,
                        value = 50),
            submitButton("Submit")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            fluidRow(
                box(
                    plotOutput("R2_prior")
                ),
                box(
                    plotOutput("w_prior")
                ),
                box(
                    plotOutput("phi_prior")
                ),
                box(
                    plotOutput("beta_prior")
                )
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    sample_R2D2 <- reactive({
        R2 <- rbeta(input$N,shape1=input$a,shape2 = input$b)
        w <- rbetapr(input$N,shape1=input$a,shape2=input$b)
        phi <- rdirichlet(input$N,alpha=rep(input$a_pi,input$p))
        sigma <- 1
        #sigma <- rexp(input$N,rate=1)
        beta <- sapply(1:input$p,function(j) rnorm(input$N,0,sd=sqrt(2)*sigma*sqrt(0.5*phi[,j]*w)))
        return(list(R2=R2,w=w,phi=phi,sigma=sigma,beta=beta))
    })
    output$R2_prior <- renderPlot({
        sample_list <- sample_R2D2()
        ggplot(data=data.frame(R2=sample_list$R2),aes(R2)) %>%
        gg_histogram_style(expression(paste(R^{2}," - prior distribution")),input$bins)
    })
    output$w_prior <- renderPlot({
        sample_list <- sample_R2D2()
        ggplot(data=data.frame(w=sample_list$w),aes(w)) %>%
            gg_histogram_style(expression(paste(omega," - prior distribution")),input$bins)
    })
    output$phi_prior <- renderPlot({
        sample_list <- sample_R2D2()
        ggplot(data=data.frame(phi=sample_list$phi[,1]),aes(phi)) %>%
            gg_histogram_style(expression(paste(phi[j]," - prior distribution")),input$bins)
    })
    output$beta_prior <- renderPlot({
        sample_list <- sample_R2D2()
        ggplot(data=data.frame(beta=sample_list$beta[,1]),aes(beta)) %>%
            gg_histogram_style(expression(paste(beta[j]," - prior distribution")),input$bins)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
