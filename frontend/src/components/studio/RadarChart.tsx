
import React from 'react';
import {
    Chart as ChartJS,
    RadialLinearScale,
    PointElement,
    LineElement,
    Filler,
    Tooltip,
    Legend,
} from 'chart.js';
import { Radar } from 'react-chartjs-2';

ChartJS.register(
    RadialLinearScale,
    PointElement,
    LineElement,
    Filler,
    Tooltip,
    Legend
);

interface RadarChartProps {
    scores: {
        drug_likeness: number;
        complexity: number;
        lipinski?: number;
        solubility: number;
        size: number;
        [key: string]: number | undefined;
    };
}

export const RadarChart: React.FC<RadarChartProps> = ({ scores }) => {
    // Normalize nomenclature for display
    const data = {
        labels: ['Drug-likeness', 'Complexity', 'Stability', 'Solubility', 'Optimal Size'],
        datasets: [
            {
                label: 'Molecular Profile',
                data: [
                    scores.drug_likeness,
                    scores.complexity,
                    scores.lipinski || scores.lipophilicity || 0,
                    scores.solubility,
                    scores.size,
                ],
                backgroundColor: 'rgba(59, 130, 246, 0.2)',
                borderColor: 'rgb(59, 130, 246)',
                borderWidth: 2,
                pointBackgroundColor: 'rgb(59, 130, 246)',
                pointBorderColor: '#fff',
                pointHoverBackgroundColor: '#fff',
                pointHoverBorderColor: 'rgb(59, 130, 246)',
            },
        ],
    };

    const options = {
        scales: {
            r: {
                angleLines: {
                    display: true,
                    color: 'rgba(0, 0, 0, 0.05)',
                },
                suggestedMin: 0,
                suggestedMax: 1,
                ticks: {
                    display: false,
                    stepSize: 0.2,
                },
                grid: {
                    color: 'rgba(0,0,0,0.05)'
                },
                pointLabels: {
                    font: {
                        size: 9,
                        weight: 'bold' as const,
                        family: 'Inter',
                    },
                    color: '#94a3b8',
                }
            },
        },
        plugins: {
            legend: {
                display: false,
            },
            tooltip: {
                enabled: true,
                backgroundColor: 'rgba(255, 255, 255, 0.9)',
                titleColor: '#1e293b',
                bodyColor: '#1e293b',
                borderColor: '#e2e8f0',
                borderWidth: 1,
                padding: 8,
                bodyFont: { size: 10 },
                displayColors: false
            }
        },
        maintainAspectRatio: false,
    };

    return (
        <div className="w-full h-40">
            <Radar data={data} options={options} />
        </div>
    );
};
