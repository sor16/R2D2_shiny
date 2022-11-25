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

geom_line_custom <- function(p,title){
    p + geom_line(color="#000000",size=1.5) +
        geom_area(fill="#0099F8") +
        scale_y_continuous(expand=c(0, 0.02)) +
        labs(title=title,x="",y="Density") +
        gg_theme()
}

geom_histogram_custom <- function(p,title,n_bins){
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
                   h5('Prior concentration'),
                   sliderInput("a_pi",
                               "a_pi",
                               value=1, 
                               min=0,max=10,step=0.01),
                   textInput('a_pi_text','Prior concentration (custom value)')
               ),
               wellPanel(
                   h4("Data setup"),
                   sliderInput("n",
                               "Sample size",
                               value=100, 
                               min = 100,max=1000,step=1),
                   sliderInput("p",
                               "Number of regression parameters",
                               value=4, 
                               min = 1,max=100,step=1),
               ),
               wellPanel(
                   h4("Plot settings"),
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
                    #plotOutput("tau2_prior")
                    plotOutput("phi_prior")
                ),
                box(
                    plotOutput("p_eff_prior")
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
        phi <- rdirichlet(input$n_samples,alpha=rep(vals$a_pi,input$p))
        sigma <- rexp(input$n_samples,rate=1)
        beta <- sapply(1:input$p,function(j) rnorm(input$n_samples,0,sd=sqrt(2)*sigma*sqrt(0.5*phi[,j]*tau2)))
        #shrinkage factor
        SS_x=matrix(rchisq(input$n_samples*input$p,df = input$n),nrow=input$n_samples)
        kappa <-  1 / (1 + tau2*phi*SS_x)
        #num effective parameters
        p_eff <- rowSums(1-kappa)
        return(list(R2=R2,tau2=tau2,phi=phi,sigma=sigma,beta=beta,p_eff=p_eff))
    })
    
    vals <- reactiveValues(a_pi = 1)
    
    muphi_input <- reactive({
        list(input$mu,input$phi,input$phi_text)
    })
    
    ab_input <- reactive({
        list(input$a,input$b)
    })
    
    a_pi_input <- reactive({
        list(input$a_pi,input$a_pi_text)
    })
    
    observeEvent(muphi_input(), {
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
    
    observeEvent(ab_input(), {
        if(input$param=='ab'){
            muphi <- get_muphi_beta_param(input$a,input$b)
            updateSliderInput(session, "phi", max=min(1000,ceiling(muphi$phi)))
            updateSliderInput(session, "mu", value=muphi$mu)
            updateSliderInput(session, "phi", value=muphi$phi)
        }
    })
    
    observeEvent(a_pi_input(), {
        if(is.na(as.numeric(input$a_pi_text))){
            a_pi_text <- ''
        }else if(as.numeric(input$a_pi_text)<=0){
            a_pi_text <- ''
        }else{
            a_pi_text <- input$a_pi_text
        }
        print(input$a_pi_text)
        print(a_pi_text)
        if(nchar(a_pi_text)==0){
            vals$a_pi <- input$a_pi
        }else{
            vals$a_pi <- as.numeric(a_pi_text)
        }
        print(vals$a_pi)
        
    })
    output$R2_prior <- renderPlot({
        #marginal dirichlet for parameter j is beta with shape1=a_j and shape2=sum_{i!=j} a_i 
        plot_dat=tibble(x=seq(0,1,by=0.001),R2=dbeta(x,shape1=input$a,shape2=input$b))
        ggplot(plot_dat,aes(x,R2)) %>%
        geom_line_custom(expression(paste(R^{2}," - prior distribution (~Beta(a,b))"))) +
        scale_x_continuous(limits=c(0,1))
    })
    output$phi_prior <- renderPlot({
        #marginal dirichlet for parameter j is beta with shape1=a_j and shape2=sum_{i!=j} a_i 
        plot_dat=tibble(x=seq(0,1,by=0.001),phi=dbeta(x,shape1=vals$a_pi,shape2=(input$p-1)*vals$a_pi))
        ggplot(data=plot_dat,aes(x,phi)) %>%
            geom_line_custom(expression(paste(phi[j]," - marginal prior distribution (",phi,"~","Dirichlet(",a[pi],"...",a[pi],"))"))) +
            scale_x_continuous(limits=c(0,1))
    })
    
    output$p_eff_prior <- renderPlot({
        sample_list <- sample_R2D2()
        ggplot(data=data.frame(p_eff_prop=sample_list$p_eff/input$p),aes(p_eff_prop)) %>%
            geom_histogram_custom(expression(paste(p[eff]/p," - prior distribution")),input$bins) +
            scale_x_continuous(limits=c(0,1))
    })
    
    output$beta_prior <- renderPlot({
        sample_list <- sample_R2D2()
        ggplot(data=data.frame(beta=sample_list$beta[,1]),aes(beta)) %>%
        geom_histogram_custom(expression(paste(beta[j]," - marginal prior distribution (~N(0,",sigma^2,tau^{2},phi[j],"))")),input$bins) +
        scale_x_continuous(limits=c(-3,3))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
