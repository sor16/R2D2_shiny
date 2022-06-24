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
        labs(title=title,x="",y="Count") +
        gg_theme()
}

get_ab_beta_param <- function(mu, tau) {
    a  <- ((1 - mu)*tau - 1 / mu) * mu ^ 2
    b <- a * (1 / mu - 1)
    return(params = list(a=a,b=b))
}

get_mutau_beta_param <- function(a, b) {
    mu <- a/(a+b)
    sigma2 <- a*b/((a+b)^2*(a+b+1))
    return(params = list(mu=mu,tau=1/sigma2))
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
    fluidRow(
        column(3,
               wellPanel(
                   radioButtons("param","R2 prior parameterization",choiceNames=c("(mean,precision)","(a,b)"),choiceValues=c('mutau','ab')),
                   conditionalPanel(condition = "input.param == 'mutau'",
                           sliderInput("mu",
                                       "Mean",
                                       value=0.5, 
                                       min=0,max=1,step=0.01),
                           sliderInput("tau",
                                       "Precision",
                                       value=50, 
                                       min=4,max=100,step=1)
                    ),
                   conditionalPanel(condition = "input.param == 'ab'",
                       sliderInput("a",
                                   "a",
                                   value=get_ab_beta_param(0.5,100)$a, 
                                   min=0,max=20,step=0.1),
                       sliderInput("b",
                                   "b",
                                   value=get_ab_beta_param(0.5,100)$b, 
                                   min=0,max=20,step=0.1)
                   ),
                   h4('Prior concentration'),
                   sliderInput("a_pi",
                               "a_pi",
                               value=1, 
                               min=0,max=10,step=0.01),
               ),
               wellPanel(
                   h4("Settings"),
                   sliderInput("p",
                               "Number of regression parameters",
                               value=4, 
                               min = 1,max=100,step=1),
                   
                   sliderInput("N",
                               "Number of samples",
                               value = 1000, 
                               min = 0,max=10000,step=10),
                   sliderInput("bins",
                               "Number of bins:",
                               min = 0,
                               max = 100,
                               value = 50)
               )
               
        ),
        column(9,
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
server <- function(input, output,session) {
    sample_R2D2 <- reactive({
        R2 <- rbeta(input$N,shape1=input$a,shape2 = input$b)
        w <- rbetapr(input$N,shape1=input$a,shape2=input$b)
        phi <- rdirichlet(input$N,alpha=rep(input$a_pi,input$p))
        sigma <- 1
        #sigma <- rexp(input$N,rate=1)
        beta <- sapply(1:input$p,function(j) rnorm(input$N,0,sd=sqrt(2)*sigma*sqrt(0.5*phi[,j]*w)))
        return(list(R2=R2,w=w,phi=phi,sigma=sigma,beta=beta))
    })
    
    observeEvent(input$mu, {
        if(input$param=='mutau'){
            updateSliderInput(session, "tau", min=ceiling(1/(input$mu*(1-input$mu))))
            shape <- get_ab_beta_param(input$mu,input$tau)
            updateSliderInput(session, "a", max=max(20,ceiling(shape$a)))
            updateSliderInput(session, "b", max=max(20,ceiling(shape$b)))
            updateSliderInput(session, "a", value=shape$a)
            updateSliderInput(session, "b", value=shape$b)
        }
        
    })
    
    observeEvent(input$tau, {
        if(input$param=='mutau'){
            shape <- get_ab_beta_param(input$mu,input$tau)
            updateSliderInput(session, "a", max=max(20,ceiling(shape$a)))
            updateSliderInput(session, "b", max=max(20,ceiling(shape$b)))
            updateSliderInput(session, "a", value=shape$a)
            updateSliderInput(session, "b", value=shape$b)
        }

    })
    
    observeEvent(input$a, {
        if(input$param=='ab'){
            mutau <- get_mutau_beta_param(input$a,input$b)
            updateSliderInput(session, "tau", max=max(1000,ceiling(mutau$tau)))
            updateSliderInput(session, "mu", value=mutau$mu)
            updateSliderInput(session, "tau", value=mutau$tau)
        }
    })
    
    observeEvent(input$b, {
        if(input$param=='ab'){
            mutau <- get_mutau_beta_param(input$a,input$b)
            updateSliderInput(session, "tau", max=max(1000,ceiling(mutau$tau)))
            updateSliderInput(session, "mu", value=mutau$mu)
            updateSliderInput(session, "tau", value=mutau$tau)
        }
    })
    output$R2_prior <- renderPlot({
        sample_list <- sample_R2D2()
        ggplot(data=data.frame(R2=sample_list$R2),aes(R2)) %>%
        gg_histogram_style(expression(paste(R^{2}," - prior distribution (~Beta(a,b))")),input$bins) +
        scale_x_continuous(limits=c(0,1))
    })
    output$w_prior <- renderPlot({
        sample_list <- sample_R2D2()
        ggplot(data=data.frame(w=sample_list$w),aes(w)) %>%
            gg_histogram_style(expression(paste(omega," - prior distribution (~BetaPrime(a,b))")),input$bins)
    })
    output$phi_prior <- renderPlot({
        sample_list <- sample_R2D2()
        ggplot(data=data.frame(phi=sample_list$phi[,1]),aes(phi)) %>%
            gg_histogram_style(expression(paste(phi[j]," - marginal prior distribution (",phi,"~","Dirichlet(",a[pi],"...",a[pi],"))")),input$bins) +
            scale_x_continuous(limits=c(0,1))
    })
    output$beta_prior <- renderPlot({
        sample_list <- sample_R2D2()
        ggplot(data=data.frame(beta=sample_list$beta[,1]),aes(beta)) %>%
        gg_histogram_style(expression(paste(beta[j]," - prior distribution (~N(0,",sigma^2,phi[j],omega,"))")),input$bins) +
        scale_x_continuous(limits=c(-3,3))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
