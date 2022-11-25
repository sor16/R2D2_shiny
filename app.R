library(shiny)
library(extraDistr)
library(ggplot2)
library(shinydashboard)
library(dplyr)
library(simstudy)

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

get_ab_beta_param <- function(mu, phi) {
    shapes <- betaGetShapes(mu,phi)
    return(params = list(a=shapes$shape1,b=shapes$shape2))
}

get_muphi_beta_param <- function(a, b) {
    mu <- a/(a+b)
    phi <- a+b
    return(params = list(mu=mu,phi=phi))
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    withMathJax(),
    # Application title
    titlePanel("R2-D2 prior distribution"),
    fluidRow(
        column(3,
               wellPanel(
                   radioButtons("param","R2 prior parameterization",choiceNames=c("(mean,precision)","(a,b)"),choiceValues=c('muphi','ab')),
                   conditionalPanel(condition = "input.param == 'muphi'",
                           sliderInput("mu",
                                       "Mean",
                                       value=0.5, 
                                       min=0,max=1,step=0.01),
                           sliderInput("phi",
                                       "Precision",
                                       value=10, 
                                       min=0.5,max=20,step=0.5),
                           textInput('phi_text','Precision (custom value)')
                    ),
                   conditionalPanel(condition = "input.param == 'ab'",
                       sliderInput("a",
                                   "a",
                                   value=get_ab_beta_param(0.5,10)$a, 
                                   min=0,max=20,step=0.1),
                       sliderInput("b",
                                   "b",
                                   value=get_ab_beta_param(0.5,10)$b, 
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
                   sliderInput("n",
                               "Sample size",
                               value=100, 
                               min = 100,max=1000,step=1),
                   sliderInput("p",
                               "Number of regression parameters",
                               value=4, 
                               min = 1,max=100,step=1),
                   
                   sliderInput("n_samples",
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
                    plotOutput("tau2_prior")
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
        set.seed(1)
        R2 <- rbeta(input$n_samples,shape1=input$a,shape2 = input$b)
        tau2 <- rbetapr(input$n_samples,shape1=input$a,shape2=input$b)
        phi <- rdirichlet(input$n_samples,alpha=rep(input$a_pi,input$p))
        #sigma <- 1
        sigma <- rexp(input$n_samples,rate=1)
        beta <- sapply(1:input$p,function(j) rnorm(input$n_samples,0,sd=sqrt(2)*sigma*sqrt(0.5*phi[,j]*tau2)))
        #shrinkage factor, sample
        SS_x=matrix(rchisq(input$n_samples*input$p,df = input$n),nrow=input$n_samples)
        kappa <-  1 / (1 + tau2*phi*SS_x)
        p_eff <- colSums(1-kappa)
        return(list(R2=R2,tau2=tau2,phi=phi,sigma=sigma,beta=beta,p_eff=p_eff))
    })
    
    muphi <- reactive({
        list(input$mu,input$phi,input$phi_text)
    })
    
    ab <- reactive({
        list(input$a,input$b)
    })
    
    observeEvent(muphi(), {
        if(input$param=='muphi'){
            if(is.na(as.numeric(input$phi_text))){
                phi_text <- ''
            }else if(as.numeric(input$phi_text)<=0){
                phi_text <- ''
            }else{
                phi_text <- input$phi_text
            }
            shape <- get_ab_beta_param(input$mu,if(nchar(phi_text)==0) input$phi else as.numeric(phi_text))
            updateSliderInput(session, "a", max=max(20,ceiling(shape$a)))
            updateSliderInput(session, "b", max=max(20,ceiling(shape$b)))
            updateSliderInput(session, "a", value=shape$a)
            updateSliderInput(session, "b", value=shape$b)
        }
        
    })
    
    observeEvent(ab(), {
        if(input$param=='ab'){
            muphi <- get_muphi_beta_param(input$a,input$b)
            updateSliderInput(session, "phi", max=min(1000,ceiling(muphi$phi)))
            updateSliderInput(session, "mu", value=muphi$mu)
            updateSliderInput(session, "phi", value=muphi$phi)
        }
    })
    output$R2_prior <- renderPlot({
        sample_list <- sample_R2D2()
        ggplot(data=data.frame(R2=sample_list$R2),aes(R2)) %>%
        gg_histogram_style(expression(paste(R^{2}," - prior distribution (~Beta(a,b))")),input$bins) +
        scale_x_continuous(limits=c(0,1))
    })
    output$tau2_prior <- renderPlot({
        sample_list <- sample_R2D2()
        ggplot(data=data.frame(tau2=sample_list$tau2),aes(tau2)) %>%
            gg_histogram_style(expression(paste(tau^{2}," - prior distribution (~BetaPrime(a,b))")),input$bins)
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
        gg_histogram_style(expression(paste(beta[j]," - marginal prior distribution (~N(0,",sigma^2,phi[j],omega,"))")),input$bins) +
        scale_x_continuous(limits=c(-3,3))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
